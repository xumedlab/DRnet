import csv
import runpy

import pandas as pd

from pipeline_utils import ensure_dirs, log_message, read_gmt

cfg = runpy.run_path('00_config.py')
RAW_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['RESULT_DIR']


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    hallmark = set(read_gmt(RAW_DIR / 'hallmark_inflammatory_response.gmt')['HALLMARK_INFLAMMATORY_RESPONSE'])

    deg = pd.read_csv(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    deg_sig = deg[(deg['significant'] == 1) & deg['gene_symbol'].fillna('').ne('')].copy()
    core = deg_sig[deg_sig['gene_symbol'].isin(hallmark)].copy()

    trend = pd.read_csv(RESULT_DIR / 'tables' / 'severity_trend_all_genes.csv')
    trend = trend[['gene', 'spearman_rho', 'padj', 'trend_significant']].rename(
        columns={'padj': 'trend_padj'}
    )
    merged = core.merge(trend, on='gene', how='left')
    progressive = merged[merged['trend_significant'] == 1].copy()

    core_out = core[
        ['gene', 'gene_symbol', 'baseMean', 'log2FC', 'pvalue', 'padj']
    ].rename(columns={'gene': 'ensembl_id'})
    progressive_out = progressive[
        ['gene', 'gene_symbol', 'log2FC', 'padj', 'spearman_rho', 'trend_padj']
    ].rename(columns={'gene': 'ensembl_id', 'padj': 'deg_padj', 'spearman_rho': 'rho'})

    core_out.to_csv(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv', index=False)
    progressive_out.to_csv(RESULT_DIR / 'tables' / 'progressive_inflammatory_genes.csv', index=False)

    summary_rows = [
        {'metric': 'primary_deg_significant', 'value': int(len(deg_sig))},
        {'metric': 'hallmark_inflammatory_genes', 'value': int(len(hallmark))},
        {'metric': 'inflammatory_core_genes', 'value': int(len(core_out))},
        {'metric': 'progressive_inflammatory_genes', 'value': int(len(progressive_out))},
    ]
    with open(RESULT_DIR / 'tables' / 'candidate_gene_summary.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['metric', 'value'])
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        '# Candidate mapping check',
        '',
        f'- significant primary DEGs with mapped symbol: {len(deg_sig)}',
        f'- hallmark inflammatory response genes: {len(hallmark)}',
        f'- inflammatory core genes: {len(core_out)}',
        f'- progressive inflammatory genes: {len(progressive_out)}',
    ]
    (RESULT_DIR / 'logs' / '06_candidate_selection_mapping_check.md').write_text(
        '\n'.join(lines) + '\n',
        encoding='utf-8',
    )

    log_message(
        '06_candidate_selection',
        f'core={len(core_out)} progressive={len(progressive_out)}',
    )


if __name__ == '__main__':
    main()
