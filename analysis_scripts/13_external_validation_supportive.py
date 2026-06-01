import csv
import gzip
import math
import runpy
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from pipeline_utils import bh_adjust, ensure_dirs, log_message

cfg = runpy.run_path('00_config.py')
RESULT_DIR = cfg['RESULT_DIR']

GSE102485_FPKM = Path('GSE102485_expressed_gene_FPKM.txt.gz')
GSE102485_META = Path('GSE102485_series_matrix.txt.gz')
GSE53257_META = Path('GSE53257_series_matrix.txt.gz')
GPL18056_META = Path('GPL18056_full.txt')


def mann_whitney_u(x, y):
    vals = [(float(v), 0) for v in x] + [(float(v), 1) for v in y]
    vals.sort(key=lambda z: z[0])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j < len(vals) and vals[j][0] == vals[i][0]:
            j += 1
        rank_value = (i + j + 1) / 2.0
        for k in range(i, j):
            ranks[k] = rank_value
        i = j
    rx = sum(r for r, (_, g) in zip(ranks, vals) if g == 0)
    nx, ny = len(x), len(y)
    u = rx - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    sigma = np.sqrt(nx * ny * (nx + ny + 1) / 12.0) if nx * ny > 0 else 1.0
    z = (u - mu) / (sigma + 1e-9)
    p = 2.0 * (1.0 - 0.5 * (1.0 + math.erf(abs(z) / np.sqrt(2.0))))
    return float(u), float(p)


def parse_series_matrix_metadata(path):
    rows = {}
    with gzip.open(path, 'rt', encoding='utf-8', errors='ignore') as handle:
        for line in handle:
            if not line.startswith('!Sample_'):
                continue
            rec = next(csv.reader([line], delimiter='\t'))
            key = rec[0].replace('!Sample_', '')
            values = [x.strip().strip('"') for x in rec[1:]]
            rows.setdefault(key, []).append(values)

    samples = rows['title'][0]
    data = [{'sample_id': sample_id} for sample_id in samples]
    for key, value_lists in rows.items():
        if key == 'title':
            continue
        if key == 'characteristics_ch1':
            for values in value_lists:
                for idx, value in enumerate(values):
                    if ': ' in value:
                        inner_key, inner_value = value.split(': ', 1)
                        data[idx][inner_key.strip()] = inner_value.strip()
            continue
        values = value_lists[0]
        for idx, value in enumerate(values):
            data[idx][key] = value
    return pd.DataFrame(data)


def load_gse102485_expression():
    df = pd.read_csv(GSE102485_FPKM, sep='\t', compression='gzip')
    sample_cols = df.columns[1:31].tolist()
    df = df.rename(columns={'Symbol': 'gene_symbol'})
    df['gene_symbol'] = df['gene_symbol'].astype(str).str.upper()
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors='coerce')
    df['mean_expr'] = df[sample_cols].mean(axis=1)
    df = (
        df[df['gene_symbol'].notna() & df['gene_symbol'].ne('--')]
        .sort_values('mean_expr', ascending=False)
        .drop_duplicates('gene_symbol')
        .set_index('gene_symbol')
    )
    return df, sample_cols


def build_gse102485_groups(meta):
    out = meta.copy()
    out['group'] = 'other'
    diabetic_mask = out['disease'].str.contains('type I diabetes|type II diabetes', case=False, regex=True)
    membrane_mask = out['source_name_ch1'].str.contains('proliferative membrane', case=False)
    normal_mask = out['disease'].str.contains('Normal retina', case=False)
    non_diabetic_disease_mask = out['disease'].str.contains(
        'branch retinal vein occlusion|retinal periphlebitis',
        case=False,
        regex=True,
    )

    out.loc[diabetic_mask & membrane_mask, 'group'] = 'diabetic_proliferative_membrane'
    out.loc[non_diabetic_disease_mask & membrane_mask, 'group'] = 'non_diabetic_membrane_control'
    out.loc[normal_mask, 'group'] = 'normal_retina'
    return out


def build_dataset_screening():
    signature = pd.read_csv(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv')['gene_symbol'].str.upper().tolist()
    gpl = pd.read_csv(GPL18056_META, sep='\t', comment='!', skiprows=44)
    gpl['Gene Symbol'] = gpl['Gene Symbol'].astype(str).str.upper()
    gse53257_overlap = sorted(set(signature).intersection(set(gpl['Gene Symbol'])))

    gse102485 = pd.read_csv(GSE102485_FPKM, sep='\t', compression='gzip')
    gse102485['Symbol'] = gse102485['Symbol'].astype(str).str.upper()
    gse102485_overlap = sorted(set(signature).intersection(set(gse102485['Symbol'])))

    rows = [
        {
            'dataset': 'GSE53257',
            'platform_or_type': 'custom mitoscriptome microarray',
            'signature_gene_overlap': len(gse53257_overlap),
            'status': 'screened_out',
            'note': 'Insufficient overlap with the current seven-gene signature; unsuitable as a supportive external-comparison dataset.',
        },
        {
            'dataset': 'GSE102485',
            'platform_or_type': 'RNA-seq of proliferative membrane / retina samples',
            'signature_gene_overlap': len(gse102485_overlap),
            'status': 'used_as_supportive_comparison',
            'note': 'Used as a supportive cross-compartment external comparison for advanced proliferative DR.',
        },
    ]
    out_path = RESULT_DIR / 'tables' / 'external_validation_dataset_screening.csv'
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return rows


def compute_signature_scores(expr_df, sample_ids, coefficients):
    available = [gene for gene in coefficients.index if gene in expr_df.index]
    matrix = np.log2(expr_df.loc[available, sample_ids].astype(float) + 1.0)
    matrix = matrix.T
    means = matrix.mean(axis=0)
    stds = matrix.std(axis=0, ddof=0).replace(0.0, 1.0)
    z = (matrix - means) / stds
    score = z.mul(coefficients.loc[available], axis=1).sum(axis=1)
    return score


def make_figure(meta, signature_df, gene_stats):
    palette = {
        'diabetic_proliferative_membrane': '#c23b4d',
        'non_diabetic_membrane_control': '#4b6cb7',
        'normal_retina': '#4b9b61',
    }
    display = {
        'diabetic_proliferative_membrane': 'Diabetic membrane',
        'non_diabetic_membrane_control': 'Non-diabetic membrane control',
        'normal_retina': 'Normal retina',
    }
    order = ['diabetic_proliferative_membrane', 'non_diabetic_membrane_control', 'normal_retina']

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    ax = axes[0]
    positions = np.arange(1, len(order) + 1)
    values = [signature_df[signature_df['group'] == group]['signature_score'].tolist() for group in order]
    box = ax.boxplot(values, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
    for patch, group in zip(box['boxes'], order):
        patch.set_facecolor(palette[group])
        patch.set_alpha(0.35)
    rng = np.random.default_rng(20250410)
    for pos, group, vals in zip(positions, order, values):
        jitter = rng.normal(0, 0.05, size=len(vals))
        ax.scatter(
            np.full(len(vals), pos) + jitter,
            vals,
            s=28,
            color=palette[group],
            alpha=0.85,
            edgecolor='white',
            linewidth=0.4,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels([display[group] for group in order], rotation=15, ha='right')
    ax.set_ylabel('External signature score')
    ax.set_title('Supportive cross-compartment comparison in GSE102485')
    ax.grid(axis='y', alpha=0.2)

    ax = axes[1]
    plot_df = gene_stats.sort_values('delta_dm_vs_nondm_mem', ascending=False).copy()
    y = np.arange(len(plot_df))
    colors = ['#c23b4d' if value >= 0 else '#2f6fb0' for value in plot_df['delta_dm_vs_nondm_mem']]
    ax.barh(y, plot_df['delta_dm_vs_nondm_mem'], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df['gene_symbol'])
    ax.set_xlabel('Delta log2(FPKM+1): diabetic membrane - non-diabetic membrane control')
    ax.set_title('Seven-gene direction in GSE102485')
    ax.axvline(0.0, color='black', linewidth=0.8)
    max_delta = float(plot_df['delta_dm_vs_nondm_mem'].max())
    ax.set_xlim(0.0, max_delta + 0.16)
    ax.grid(axis='x', alpha=0.2)
    for yi, (_, row) in zip(y, plot_df.iterrows()):
        ax.text(
            row['delta_dm_vs_nondm_mem'] + (0.03 if row['delta_dm_vs_nondm_mem'] >= 0 else -0.03),
            yi,
            f"FDR={row['padj_dm_vs_nondm_mem']:.3f}",
            va='center',
            ha='left' if row['delta_dm_vs_nondm_mem'] >= 0 else 'right',
            fontsize=7,
        )

    fig.tight_layout()
    for ext in cfg['FIG_FORMATS']:
        fig.savefig(RESULT_DIR / 'figures' / f'Frontiers_supportive_external_validation.{ext}', dpi=cfg['FIG_DPI'], bbox_inches='tight')
    plt.close(fig)


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'figures', RESULT_DIR / 'logs')
    screening_rows = build_dataset_screening()

    expr_df, sample_cols = load_gse102485_expression()
    meta = build_gse102485_groups(parse_series_matrix_metadata(GSE102485_META))
    meta_out = meta[['sample_id', 'group', 'disease', 'source_name_ch1', 'tissue']].copy()
    meta_out.to_csv(RESULT_DIR / 'tables' / 'external_validation_gse102485_metadata.csv', index=False)

    coefficients = pd.read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    coefficients['gene_symbol'] = coefficients['gene_symbol'].str.upper()
    coef_series = coefficients.set_index('gene_symbol')['coefficient']

    signature_scores = compute_signature_scores(expr_df, sample_cols, coef_series)
    signature_df = meta[['sample_id', 'group']].merge(
        signature_scores.rename('signature_score'),
        left_on='sample_id',
        right_index=True,
        how='left',
    )
    diabetic_scores = signature_df[signature_df['group'] == 'diabetic_proliferative_membrane']['signature_score']
    non_dm_scores = signature_df[signature_df['group'] == 'non_diabetic_membrane_control']['signature_score']
    if len(set(np.concatenate([np.ones(len(diabetic_scores)), np.zeros(len(non_dm_scores))]))) >= 2:
        y_true = np.array([1] * len(diabetic_scores) + [0] * len(non_dm_scores))
        y_score = np.concatenate([diabetic_scores.to_numpy(dtype=float), non_dm_scores.to_numpy(dtype=float)])
        signature_auc = float(roc_auc_score(y_true, y_score))
    else:
        signature_auc = float('nan')
    signature_df.to_csv(RESULT_DIR / 'tables' / 'external_validation_gse102485_signature_scores.csv', index=False)

    selected_genes = coefficients['gene_symbol'].tolist()
    rows = []
    pvals = []
    for gene_symbol in selected_genes:
        if gene_symbol not in expr_df.index:
            continue
        values = np.log2(expr_df.loc[gene_symbol, sample_cols].astype(float) + 1.0)
        gene_df = pd.DataFrame({'sample_id': sample_cols, 'expression': values.to_numpy(dtype=float)})
        gene_df = gene_df.merge(meta[['sample_id', 'group']], on='sample_id', how='left')
        diabetic_mem = gene_df[gene_df['group'] == 'diabetic_proliferative_membrane']['expression']
        non_dm_mem = gene_df[gene_df['group'] == 'non_diabetic_membrane_control']['expression']
        normal_retina = gene_df[gene_df['group'] == 'normal_retina']['expression']
        _, pvalue = mann_whitney_u(diabetic_mem.tolist(), non_dm_mem.tolist())
        pvals.append(pvalue)
        rows.append(
            {
                'gene_symbol': gene_symbol,
                'mean_dm_mem': float(diabetic_mem.mean()),
                'mean_non_dm_mem': float(non_dm_mem.mean()),
                'mean_normal_retina': float(normal_retina.mean()),
                'delta_dm_vs_nondm_mem': float(diabetic_mem.mean() - non_dm_mem.mean()),
                'delta_dm_vs_normal_retina': float(diabetic_mem.mean() - normal_retina.mean()),
                'pvalue_dm_vs_nondm_mem': float(pvalue),
                'direction_match_primary': int((diabetic_mem.mean() - non_dm_mem.mean()) > 0),
            }
        )
    padj = bh_adjust(pvals)
    for row, adj in zip(rows, padj):
        row['padj_dm_vs_nondm_mem'] = float(adj)
    gene_stats = pd.DataFrame(rows).sort_values(['padj_dm_vs_nondm_mem', 'delta_dm_vs_nondm_mem'], ascending=[True, False])
    gene_stats.to_csv(RESULT_DIR / 'tables' / 'external_validation_gse102485_gene_stats.csv', index=False)

    summary = pd.DataFrame(
        [
            {
                'metric': 'screened_datasets',
                'value': len(screening_rows),
            },
            {
                'metric': 'gse102485_diabetic_proliferative_membrane_n',
                'value': int((meta['group'] == 'diabetic_proliferative_membrane').sum()),
            },
            {
                'metric': 'gse102485_non_diabetic_membrane_control_n',
                'value': int((meta['group'] == 'non_diabetic_membrane_control').sum()),
            },
            {
                'metric': 'gse102485_normal_retina_n',
                'value': int((meta['group'] == 'normal_retina').sum()),
            },
            {
                'metric': 'supportive_signature_auc_dm_mem_vs_nondm_mem',
                'value': signature_auc,
            },
            {
                'metric': 'genes_with_positive_direction_dm_mem_vs_nondm_mem',
                'value': int((gene_stats['delta_dm_vs_nondm_mem'] > 0).sum()),
            },
        ]
    )
    summary.to_csv(RESULT_DIR / 'tables' / 'external_validation_supportive_summary.csv', index=False)

    make_figure(meta, signature_df, gene_stats)
    log_message(
        '13_external_validation_supportive',
        (
            f"screened_out=GSE53257 used=GSE102485 "
            f"signature_auc={signature_auc:.4f} positive_direction={(gene_stats['delta_dm_vs_nondm_mem'] > 0).sum()}"
        ),
    )


if __name__ == '__main__':
    main()
