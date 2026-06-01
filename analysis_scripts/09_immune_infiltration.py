import csv
import runpy

import gseapy as gp
import pandas as pd

from pipeline_utils import (
    bh_adjust,
    ensure_dirs,
    load_ensembl_symbol_mapping,
    log_message,
    mann_whitney_u,
    matrix_ensembl_to_symbol,
    read_gmt,
    spearman,
)

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
SEED = cfg.get('RANDOM_SEED', 202501)


def build_symbol_expression():
    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        matrix_ensembl = {row[0]: [float(x) for x in row[1:]] for row in reader}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    matrix_symbol, _ = matrix_ensembl_to_symbol(matrix_ensembl, mapping)
    expr = pd.DataFrame(matrix_symbol, index=header[1:]).T
    expr.index.name = 'gene_symbol'
    return header[1:], expr


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    gene_sets = read_gmt(PROC_DIR / 'immune_28_signatures_charoentong2017_human_cleaned.gmt')
    samples, expr = build_symbol_expression()

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        pheno = {row['sample_id']: row for row in csv.DictReader(f)}

    ss = gp.ssgsea(
        data=expr,
        gene_sets=gene_sets,
        outdir=None,
        sample_norm_method='rank',
        correl_norm_type='rank',
        min_size=3,
        max_size=5000,
        weight=0.25,
        threads=1,
        no_plot=True,
        seed=SEED,
        verbose=False,
    )
    score_df = ss.res2d.rename(columns={'Name': 'sample_id', 'Term': 'cell_type', 'NES': 'score'})[
        ['sample_id', 'cell_type', 'score']
    ].copy()
    score_df['score'] = score_df['score'].astype(float)
    score_df['group'] = score_df['sample_id'].map(lambda x: pheno[x]['disease_group'])
    score_df['severity_code'] = score_df['sample_id'].map(lambda x: pheno[x]['severity_code'])
    score_df = score_df.sort_values(['cell_type', 'sample_id'])
    score_df.to_csv(RESULT_DIR / 'tables' / 'immune_ssgsea_scores.csv', index=False)

    comp_rows = []
    trend_rows = []
    for cell_type, sub_df in score_df.groupby('cell_type'):
        ctrl = sub_df[sub_df['group'] == cfg['PRIMARY_CTRL']]['score'].tolist()
        case = sub_df[sub_df['group'] == cfg['PRIMARY_CASE']]['score'].tolist()
        _, pvalue = mann_whitney_u(ctrl, case)
        comp_rows.append(
            {
                'cell_type': cell_type,
                'control_mean': sum(ctrl) / len(ctrl),
                'case_mean': sum(case) / len(case),
                'delta_mean': (sum(case) / len(case)) - (sum(ctrl) / len(ctrl)),
                'pvalue': pvalue,
            }
        )
        rho, trend_p = spearman(
            sub_df['severity_code'].astype(int).tolist(),
            sub_df['score'].tolist(),
        )
        trend_rows.append({'cell_type': cell_type, 'rho': rho, 'pvalue': trend_p})

    for rows in [comp_rows, trend_rows]:
        padj = bh_adjust([row['pvalue'] for row in rows])
        for row, adj in zip(rows, padj):
            row['padj'] = adj
            row['significant'] = int(adj < 0.05)

    pd.DataFrame(comp_rows).sort_values(['padj', 'pvalue']).to_csv(
        RESULT_DIR / 'tables' / 'immune_primary_comparison.csv',
        index=False,
    )
    pd.DataFrame(trend_rows).sort_values(['padj', 'pvalue']).to_csv(
        RESULT_DIR / 'tables' / 'immune_severity_trend.csv',
        index=False,
    )

    selected_genes = pd.read_csv(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv')['gene_symbol'].dropna().tolist()
    cor_rows = []
    for gene_symbol in selected_genes:
        if gene_symbol not in expr.index:
            continue
        gene_values = expr.loc[gene_symbol, samples].tolist()
        for cell_type, sub_df in score_df.groupby('cell_type'):
            aligned_scores = sub_df.set_index('sample_id').reindex(samples)['score'].tolist()
            rho, pvalue = spearman(gene_values, aligned_scores)
            cor_rows.append(
                {
                    'gene_symbol': gene_symbol,
                    'cell_type': cell_type,
                    'rho': rho,
                    'pvalue': pvalue,
                }
            )

    if cor_rows:
        padj = bh_adjust([row['pvalue'] for row in cor_rows])
        for row, adj in zip(cor_rows, padj):
            row['padj'] = adj
            row['significant'] = int(adj < 0.05)
    pd.DataFrame(cor_rows).sort_values(['padj', 'pvalue', 'rho'], ascending=[True, True, False]).to_csv(
        RESULT_DIR / 'tables' / 'gene_immune_correlations.csv',
        index=False,
    )

    log_message(
        '09_immune_infiltration',
        f'cell_types={len(gene_sets)} correlation_pairs={len(cor_rows)}',
    )


if __name__ == '__main__':
    main()
