import csv
import runpy

import gseapy as gp
import pandas as pd

from pipeline_utils import (
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


def load_symbol_expression():
    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        matrix_ensembl = {row[0]: [float(x) for x in row[1:]] for row in reader}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    matrix_symbol, dedup_log = matrix_ensembl_to_symbol(matrix_ensembl, mapping)
    expr = pd.DataFrame(matrix_symbol, index=header[1:]).T
    expr.index.name = 'gene_symbol'
    return header[1:], expr, dedup_log


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    hallmark_genes = set(read_gmt(RAW_DIR / 'hallmark_inflammatory_response.gmt')['HALLMARK_INFLAMMATORY_RESPONSE'])

    samples, expr, dedup_log = load_symbol_expression()
    overlap = sorted(hallmark_genes.intersection(expr.index))
    if not overlap:
        raise ValueError('No overlap between mapped gene symbols and HALLMARK_INFLAMMATORY_RESPONSE.')

    with open(PROC_DIR / 'pheno_macula_4groups.csv', encoding='utf-8') as f:
        pheno = {row['sample_id']: row for row in csv.DictReader(f)}

    ss = gp.ssgsea(
        data=expr,
        gene_sets={'HALLMARK_INFLAMMATORY_RESPONSE': overlap},
        outdir=None,
        sample_norm_method='rank',
        correl_norm_type='rank',
        min_size=5,
        max_size=5000,
        weight=0.25,
        threads=1,
        no_plot=True,
        seed=SEED,
        verbose=False,
    )
    score_df = ss.res2d.copy()
    score_df = score_df[score_df['Term'] == 'HALLMARK_INFLAMMATORY_RESPONSE'].copy()
    score_map = {row['Name']: float(row['NES']) for _, row in score_df.iterrows()}

    rows = []
    for sample in samples:
        rows.append(
            {
                'sample_id': sample,
                'group': pheno[sample]['disease_group'],
                'severity_code': pheno[sample]['severity_code'],
                'inflammation_ssgsea': score_map[sample],
            }
        )

    with open(RESULT_DIR / 'tables' / 'inflammation_ssgsea_scores.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    ctrl = [row['inflammation_ssgsea'] for row in rows if row['group'] == cfg['PRIMARY_CTRL']]
    case = [row['inflammation_ssgsea'] for row in rows if row['group'] == cfg['PRIMARY_CASE']]
    _, pvalue = mann_whitney_u(ctrl, case)
    with open(RESULT_DIR / 'tables' / 'inflammation_group_comparison.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=['contrast', 'control_mean', 'case_mean', 'delta_mean', 'pvalue'],
        )
        writer.writeheader()
        writer.writerow(
            {
                'contrast': 'healthy control vs NPDR/PDR + DME',
                'control_mean': sum(ctrl) / len(ctrl),
                'case_mean': sum(case) / len(case),
                'delta_mean': (sum(case) / len(case)) - (sum(ctrl) / len(ctrl)),
                'pvalue': pvalue,
            }
        )

    rho, trend_p = spearman(
        [int(row['severity_code']) for row in rows],
        [row['inflammation_ssgsea'] for row in rows],
    )
    with open(RESULT_DIR / 'tables' / 'inflammation_severity_trend.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['rho', 'pvalue'])
        writer.writeheader()
        writer.writerow({'rho': rho, 'pvalue': trend_p})

    mapping_log = '\n'.join(
        [
            '# Inflammation scoring mapping check',
            '',
            f'- mapped gene symbols: {expr.shape[0]}',
            f'- hallmark genes in MSigDB inflammation set: {len(hallmark_genes)}',
            f'- overlapping genes used for ssGSEA: {len(overlap)}',
            f'- deduplicated symbols kept by higher mean expression: {len(dedup_log)}',
        ]
    )
    (RESULT_DIR / 'logs' / '05_inflammation_scoring_mapping_check.md').write_text(
        mapping_log + '\n',
        encoding='utf-8',
    )
    log_message('05_inflammation_scoring', f'overlap_genes={len(overlap)}')


if __name__ == '__main__':
    main()
