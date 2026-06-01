import csv
import runpy

import gseapy as gp
import numpy as np
import pandas as pd

from pipeline_utils import (
    ensure_dirs,
    load_ensembl_symbol_mapping,
    log_message,
    matrix_ensembl_to_symbol,
    read_gmt,
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
    return expr


def run_gene_prerank(expr, gene_symbol, hallmark_sets):
    target = expr.loc[gene_symbol]
    corr = expr.apply(lambda row: row.corr(target, method='spearman'), axis=1)
    corr = corr.fillna(0.0).sort_values(ascending=False)
    jitter = np.linspace(0.0, 1e-9, len(corr), dtype=float)
    rank_df = pd.DataFrame({'gene': corr.index.astype(str), 'score': corr.values.astype(float) + jitter})
    prerank = gp.prerank(
        rnk=rank_df,
        gene_sets=hallmark_sets,
        threads=1,
        permutation_num=1000,
        min_size=5,
        max_size=5000,
        seed=SEED,
        outdir=None,
        verbose=False,
    )
    res = prerank.res2d.reset_index().rename(columns={'Term': 'pathway'})
    rows = []
    for _, row in res.iterrows():
        pvalue = float(row.get('NOM p-val', 1.0))
        padj = float(row.get('FDR q-val', 1.0))
        rows.append(
            {
                'gene_symbol': gene_symbol,
                'pathway': str(row['pathway']),
                'ES': float(row.get('ES', 0.0)),
                'NES': float(row.get('NES', 0.0)),
                'pvalue': pvalue,
                'padj': padj,
                'significant': int((pvalue < 0.05) and (padj < 0.25)),
                'geneset_size': len(hallmark_sets.get(str(row['pathway']), [])),
            }
        )
    return rows


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    selected_genes = pd.read_csv(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv')['gene_symbol'].dropna().tolist()
    if not selected_genes:
        raise ValueError('No LASSO-selected genes available for enrichment analysis.')

    expr = build_symbol_expression()
    hallmark_sets = read_gmt(RAW_DIR / 'hallmark_human_gene_symbols.gmt')

    summary_rows = []
    for gene_symbol in selected_genes:
        if gene_symbol not in expr.index:
            continue
        gene_rows = run_gene_prerank(expr, gene_symbol, hallmark_sets)
        summary_rows.extend(gene_rows)
        pd.DataFrame(gene_rows).sort_values(['padj', 'NES'], ascending=[True, False]).to_csv(
            RESULT_DIR / 'tables' / f'gsea_{gene_symbol}.csv',
            index=False,
        )

    if not summary_rows:
        raise ValueError('No GSEA rows were produced for the selected genes.')

    summary_df = pd.DataFrame(summary_rows).sort_values(['padj', 'NES'], ascending=[True, False])
    summary_df.to_csv(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv', index=False)

    sig_df = summary_df[summary_df['significant'] == 1].copy()
    if sig_df.empty:
        sig_df = summary_df.groupby('gene_symbol', as_index=False).head(5).copy()

    recurrence = (
        sig_df.groupby('pathway')
        .agg(
            n_genes=('gene_symbol', 'nunique'),
            genes=('gene_symbol', lambda x: ','.join(sorted(set(x)))),
            mean_nes=('NES', 'mean'),
            max_abs_nes=('NES', lambda x: float(abs(x).max())),
            min_padj=('padj', 'min'),
        )
        .reset_index()
        .sort_values(['n_genes', 'min_padj', 'max_abs_nes'], ascending=[False, True, False])
    )
    recurrence.to_csv(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv', index=False)

    deg = pd.read_csv(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    deg_meta = deg[['gene_symbol', 'log2FC', 'padj']].dropna().drop_duplicates('gene_symbol')
    coef = pd.read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    ipa = coef.merge(deg_meta, on='gene_symbol', how='left').sort_values('coefficient', ascending=False)
    ipa.to_csv(RESULT_DIR / 'tables' / 'ipa_input_selected_genes.csv', index=False)

    log_message(
        '08_enrichment_analysis',
        f'selected_genes={len(selected_genes)} gsea_rows={len(summary_df)} recurrent_pathways={len(recurrence)}',
    )


if __name__ == '__main__':
    main()
