import csv
import math
import runpy

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform

from pipeline_utils import ensure_dirs, load_ensembl_symbol_mapping, log_message, matrix_ensembl_to_symbol

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
FIG_DIR = RESULT_DIR / 'figures'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

GROUP_ORDER = cfg['GROUPS']
GROUP_COLORS = {
    'healthy control': '#2f6fb0',
    'diabetic': '#4b9b61',
    'NPDR': '#d98a22',
    'NPDR/PDR + DME': '#c23b4d',
}
SEED = cfg.get('RANDOM_SEED', 202501)


def display_cell_label(name):
    replacements = {
        'MDSC': 'MDSC',
    }
    if name in replacements:
        return replacements[name]
    return str(name).replace('_', ' ')


def save_figure(fig, names):
    if isinstance(names, str):
        names = [names]
    for name in names:
        for ext in cfg['FIG_FORMATS']:
            fig.savefig(FIG_DIR / f'{name}.{ext}', dpi=cfg['FIG_DPI'], bbox_inches='tight')
    plt.close(fig)


def read_csv(path):
    return pd.read_csv(path)


def load_symbol_expression():
    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        matrix_ensembl = {row[0]: [float(x) for x in row[1:]] for row in reader}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    matrix_symbol, _ = matrix_ensembl_to_symbol(matrix_ensembl, mapping)
    expr = pd.DataFrame(matrix_symbol, index=header[1:]).astype(float)
    expr.index.name = 'sample_id'
    return expr


def get_sample_order(pheno):
    ordered = []
    for group in GROUP_ORDER:
        subset = pheno[pheno['disease_group'] == group].sort_index()
        ordered.extend(subset.index.tolist())
    return ordered


def draw_boxplot(ax, df, value_col, title, ylabel):
    positions = np.arange(1, len(GROUP_ORDER) + 1)
    rng = np.random.default_rng(SEED)
    values = [df[df['group'] == group][value_col].astype(float).tolist() for group in GROUP_ORDER]
    bp = ax.boxplot(values, positions=positions, widths=0.55, patch_artist=True, showfliers=False)
    for patch, group in zip(bp['boxes'], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[group])
        patch.set_alpha(0.35)
    for pos, group, vals in zip(positions, GROUP_ORDER, values):
        jitter = rng.normal(0, 0.045, size=len(vals))
        ax.scatter(
            np.full(len(vals), pos) + jitter,
            vals,
            s=18,
            color=GROUP_COLORS[group],
            alpha=0.85,
            edgecolor='white',
            linewidth=0.3,
        )
    ax.set_xticks(positions)
    ax.set_xticklabels(['Healthy', 'Diabetic', 'NPDR', 'Advanced DR'], rotation=15)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(axis='y', alpha=0.2)


def zscore_rows(frame):
    values = frame.to_numpy(dtype=float)
    mean = values.mean(axis=1, keepdims=True)
    std = values.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    scaled = (values - mean) / std
    return pd.DataFrame(scaled, index=frame.index, columns=frame.columns)


def draw_heatmap(ax, frame, title, cmap='RdBu_r', center=0.0, show_xlabels=False):
    data = frame.to_numpy(dtype=float)
    vmax = np.nanmax(np.abs(data)) if center == 0 else np.nanmax(data)
    vmin = -vmax if center == 0 else np.nanmin(data)
    im = ax.imshow(data, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_yticks(np.arange(len(frame.index)))
    ax.set_yticklabels(frame.index, fontsize=8)
    if show_xlabels:
        ax.set_xticks(np.arange(len(frame.columns)))
        ax.set_xticklabels(frame.columns, rotation=90, fontsize=6)
    else:
        ax.set_xticks([])
    return im


def build_group_boundaries(sample_ids, pheno):
    sizes = []
    for group in GROUP_ORDER:
        sizes.append(sum(pheno.loc[sample_ids, 'disease_group'] == group))
    boundaries = np.cumsum(sizes)[:-1]
    return boundaries, sizes


def add_group_separators(ax, boundaries):
    for boundary in boundaries:
        ax.axvline(boundary - 0.5, color='black', linewidth=0.6)


def plot_workflow():
    fig, ax = plt.subplots(figsize=(12, 3.5))
    ax.axis('off')
    labels = [
        'GSE160306 macula samples',
        'Count processing\nand log2CPM',
        'DESeq2 + trend',
        'Inflammation core genes',
        'LASSO signature',
        'GSEA + immune ssGSEA',
    ]
    xs = np.linspace(0.04, 0.84, len(labels))
    for idx, (xpos, label) in enumerate(zip(xs, labels)):
        box = patches.FancyBboxPatch(
            (xpos, 0.35),
            0.12,
            0.28,
            boxstyle='round,pad=0.02,rounding_size=0.03',
            linewidth=1.2,
            edgecolor='#3a3a3a',
            facecolor='#e9eef5' if idx % 2 == 0 else '#f3efe3',
        )
        ax.add_patch(box)
        ax.text(xpos + 0.06, 0.49, label, ha='center', va='center', fontsize=10)
        if idx < len(labels) - 1:
            ax.annotate('', xy=(xs[idx + 1] - 0.01, 0.49), xytext=(xpos + 0.12, 0.49), arrowprops=dict(arrowstyle='->', lw=1.5))
    save_figure(fig, 'Figure1_workflow')


def plot_sample_composition(pheno):
    counts = pheno['disease_group'].value_counts().reindex(GROUP_ORDER)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bars = ax.bar(range(len(GROUP_ORDER)), counts.values, color=[GROUP_COLORS[g] for g in GROUP_ORDER], alpha=0.85)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(['Healthy', 'Diabetic', 'NPDR', 'Advanced DR'], rotation=15)
    ax.set_ylabel('Sample count')
    ax.set_title('Macula sample composition')
    for bar, value in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.1, str(int(value)), ha='center', va='bottom', fontsize=9)
    ax.grid(axis='y', alpha=0.2)
    save_figure(fig, 'Figure1_sample_composition')


def plot_pca():
    pca = read_csv(RESULT_DIR / 'tables' / 'pca_coordinates.csv')
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    for group in GROUP_ORDER:
        sub = pca[pca['group'] == group]
        ax.scatter(sub['PC1'], sub['PC2'], s=48, label=group, color=GROUP_COLORS[group], alpha=0.9)
    pc1 = float(pca['PC1_var_ratio'].iloc[0]) * 100.0
    pc2 = float(pca['PC2_var_ratio'].iloc[0]) * 100.0
    ax.set_xlabel(f'PC1 ({pc1:.1f}%)')
    ax.set_ylabel(f'PC2 ({pc2:.1f}%)')
    ax.set_title('PCA of macula samples')
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    save_figure(fig, ['Figure1_pca_macula_4groups', 'pca_4groups'])


def plot_inflammation():
    infl = read_csv(RESULT_DIR / 'tables' / 'inflammation_ssgsea_scores.csv')
    fig, ax = plt.subplots(figsize=(7, 4.5))
    draw_boxplot(ax, infl, 'inflammation_ssgsea', 'Inflammation ssGSEA across groups', 'Inflammation ssGSEA NES')
    save_figure(fig, 'Figure2_inflammation_ssgsea_4groups')


def plot_volcano():
    deg = read_csv(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    deg['neglog10_p'] = -np.log10(np.clip(deg['pvalue'].astype(float), 1e-300, None))
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    sig = deg['significant'] == 1
    ax.scatter(deg.loc[~sig, 'log2FC'], deg.loc[~sig, 'neglog10_p'], s=8, color='#bdbdbd', alpha=0.45)
    ax.scatter(deg.loc[sig, 'log2FC'], deg.loc[sig, 'neglog10_p'], s=10, color='#c23b4d', alpha=0.75)
    ax.axvline(cfg['LOG2FC_THRESHOLD'], linestyle='--', color='black', linewidth=0.8)
    ax.axvline(-cfg['LOG2FC_THRESHOLD'], linestyle='--', color='black', linewidth=0.8)
    ax.axhline(-math.log10(0.05), linestyle='--', color='black', linewidth=0.8)
    top = deg[deg['significant'] == 1].sort_values('padj').head(10)
    for _, row in top.iterrows():
        label = row['gene_symbol'] if isinstance(row['gene_symbol'], str) and row['gene_symbol'] else row['gene']
        ax.text(row['log2FC'], row['neglog10_p'] + 0.15, label, fontsize=7)
    ax.set_xlabel('log2 fold change')
    ax.set_ylabel('-log10(P value)')
    ax.set_title('Primary contrast volcano plot')
    ax.grid(alpha=0.15)
    save_figure(fig, 'Figure2_primary_volcano')


def plot_overlap(core_df, deg_df):
    primary_deg = int((deg_df['significant'] == 1).sum())
    hallmark_n = 200
    overlap_n = len(core_df)
    left_only = max(primary_deg - overlap_n, 0)
    right_only = max(hallmark_n - overlap_n, 0)

    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ax.axis('off')
    ax.add_patch(patches.Circle((0.42, 0.5), 0.23, color='#6d9edb', alpha=0.45))
    ax.add_patch(patches.Circle((0.58, 0.5), 0.23, color='#d98f8f', alpha=0.45))
    ax.text(0.28, 0.78, 'Primary\nDEGs', ha='center', va='center', fontsize=10)
    ax.text(0.72, 0.78, 'Hallmark\ninflammation', ha='center', va='center', fontsize=10)
    ax.text(0.33, 0.5, str(left_only), ha='center', va='center', fontsize=15, fontweight='bold')
    ax.text(0.5, 0.5, str(overlap_n), ha='center', va='center', fontsize=15, fontweight='bold')
    ax.text(0.67, 0.5, str(right_only), ha='center', va='center', fontsize=15, fontweight='bold')
    save_figure(fig, 'Figure2_overlap_primarydeg_inflammation')


def plot_gene_heatmap(expr, pheno, genes, name, title):
    genes = [gene for gene in genes if gene in expr.columns]
    if not genes:
        return
    sample_order = get_sample_order(pheno)
    matrix = expr.loc[sample_order, genes].T
    matrix = zscore_rows(matrix)
    fig, ax = plt.subplots(figsize=(10, max(4.5, 0.35 * len(genes) + 1.5)))
    im = draw_heatmap(ax, matrix, title)
    boundaries, _ = build_group_boundaries(sample_order, pheno)
    add_group_separators(ax, boundaries)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    save_figure(fig, name)


def plot_lasso_path_and_cv():
    path_df = read_csv(RESULT_DIR / 'tables' / 'lasso_path_coefficients.csv')
    coef_df = read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    keep_genes = coef_df['gene_symbol'].tolist()
    fig, ax = plt.subplots(figsize=(6.6, 5))
    for gene_symbol in keep_genes:
        sub = path_df[path_df['gene_symbol'] == gene_symbol].sort_values('lambda')
        ax.plot(np.log10(sub['lambda']), sub['coefficient'], linewidth=1.8, label=gene_symbol)
    ax.set_xlabel('log10(lambda)')
    ax.set_ylabel('Coefficient')
    ax.set_title('LASSO coefficient path')
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    save_figure(fig, 'Figure3_lasso_path')

    cv_df = read_csv(RESULT_DIR / 'tables' / 'lasso_cv_curve.csv').sort_values('lambda')
    fig, ax = plt.subplots(figsize=(6.4, 5))
    x = np.log10(cv_df['lambda'].to_numpy(dtype=float))
    y = cv_df['mean_auc'].to_numpy(dtype=float)
    yerr = cv_df['std_auc'].to_numpy(dtype=float)
    ax.plot(x, y, color='#2f6fb0', linewidth=2.0)
    ax.fill_between(x, y - yerr, y + yerr, color='#2f6fb0', alpha=0.2)
    best_idx = int(np.argmax(y))
    ax.scatter([x[best_idx]], [y[best_idx]], color='#c23b4d', s=35, zorder=3)
    ax.set_xlabel('log10(lambda)')
    ax.set_ylabel('Mean CV AUC')
    ax.set_title('Repeated CV performance')
    ax.grid(alpha=0.2)
    save_figure(fig, 'Figure3_cv_curve')


def plot_selected_gene_expression(expr, pheno):
    coef_df = read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    genes = coef_df['gene_symbol'].head(4).tolist()
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=False)
    axes = axes.flatten()
    for ax, gene_symbol in zip(axes, genes):
        sub = pd.DataFrame(
            {
                'sample_id': expr.index,
                'group': pheno.loc[expr.index, 'disease_group'].values,
                'expression': expr[gene_symbol].values,
            }
        )
        draw_boxplot(ax, sub, 'expression', gene_symbol, 'log2CPM')
    for ax in axes[len(genes):]:
        ax.axis('off')
    save_figure(fig, 'Figure3_selected_gene_expression_4groups')


def plot_roc():
    ind_curve = read_csv(RESULT_DIR / 'tables' / 'roc_curve_individual_genes.csv')
    ind_auc = read_csv(RESULT_DIR / 'tables' / 'roc_individual_genes.csv')
    fig, ax = plt.subplots(figsize=(6.3, 5.2))
    for _, row in ind_auc.sort_values('auc', ascending=False).iterrows():
        item = row['item']
        sub = ind_curve[ind_curve['item'] == item]
        ax.plot(sub['fpr'], sub['tpr'], linewidth=1.8, label=f"{item} (AUC={row['auc']:.2f})")
    ax.plot([0, 1], [0, 1], linestyle='--', color='grey', linewidth=1.0)
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('Individual gene ROC curves')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.grid(alpha=0.2)
    save_figure(fig, 'Figure3_individual_roc')

    comb_curve = read_csv(RESULT_DIR / 'tables' / 'roc_curve_combined_signature.csv')
    comb_auc = read_csv(RESULT_DIR / 'tables' / 'roc_combined_signature.csv').iloc[0]
    fig, ax = plt.subplots(figsize=(5.8, 5.0))
    ax.plot(
        comb_curve['fpr'],
        comb_curve['tpr'],
        color='#c23b4d',
        linewidth=2.2,
        label=f"Combined signature (AUC={comb_auc['auc']:.2f})",
    )
    ax.plot([0, 1], [0, 1], linestyle='--', color='grey', linewidth=1.0)
    ax.set_xlabel('False positive rate')
    ax.set_ylabel('True positive rate')
    ax.set_title('OOF ROC of combined signature')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    ax.grid(alpha=0.2)
    save_figure(fig, 'Figure3_combined_roc')


def plot_gsea():
    gsea = read_csv(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv')
    recurrence = read_csv(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv')
    top_pathways = recurrence.head(15)['pathway'].tolist()

    heat = gsea[gsea['pathway'].isin(top_pathways)].pivot(index='pathway', columns='gene_symbol', values='NES').fillna(0.0)
    heat = heat.reindex(top_pathways)
    fig, ax = plt.subplots(figsize=(7.5, max(5, 0.38 * len(heat))))
    im = draw_heatmap(ax, heat, 'Gene-specific preranked GSEA NES')
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns, rotation=45, ha='right')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    save_figure(fig, 'Figure4_selected_gene_gsea_heatmap')

    dot = recurrence.head(15).sort_values(['n_genes', 'min_padj'], ascending=[True, False])
    fdr_floor = 1e-16
    neglog_fdr_capped = -np.log10(np.clip(dot['min_padj'].astype(float), fdr_floor, 1.0))
    fig, ax = plt.subplots(figsize=(7.2, max(5, 0.35 * len(dot))))
    scatter = ax.scatter(
        dot['n_genes'],
        dot['pathway'],
        s=120 + dot['max_abs_nes'] * 55,
        c=neglog_fdr_capped,
        cmap='YlOrRd',
        vmin=0,
        vmax=-np.log10(fdr_floor),
        edgecolor='black',
        linewidth=0.4,
    )
    ax.set_xlabel('Number of signature genes with significant enrichment')
    ax.set_ylabel('Pathway')
    ax.set_title('Recurrent hallmark pathways')
    x_min = int(np.floor(dot['n_genes'].min()))
    x_max = int(np.ceil(dot['n_genes'].max()))
    ax.set_xlim(x_min - 0.45, x_max + 0.45)
    ax.set_xticks(np.arange(x_min, x_max + 1))
    fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.02, label='-log10(min FDR), capped')
    save_figure(fig, 'Figure4_recurrent_pathway_dotplot')

    selected_genes = gsea['gene_symbol'].drop_duplicates().tolist()
    net_edges = gsea[gsea['pathway'].isin(top_pathways) & gsea['gene_symbol'].isin(selected_genes)].copy()
    net_edges = net_edges.sort_values('padj').groupby(['gene_symbol', 'pathway'], as_index=False).first()
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    ax.axis('off')
    left_y = np.linspace(0.85, 0.15, len(selected_genes))
    right_y = np.linspace(0.92, 0.08, len(top_pathways))
    gene_pos = {gene: (0.18, y) for gene, y in zip(selected_genes, left_y)}
    path_pos = {pathway: (0.82, y) for pathway, y in zip(top_pathways, right_y)}
    for gene, (xpos, ypos) in gene_pos.items():
        ax.text(xpos, ypos, gene, ha='center', va='center', fontsize=10, bbox=dict(boxstyle='round,pad=0.25', fc='#dce8f5', ec='#56789a'))
    for pathway, (xpos, ypos) in path_pos.items():
        ax.text(xpos, ypos, pathway.replace('HALLMARK_', ''), ha='center', va='center', fontsize=8, bbox=dict(boxstyle='round,pad=0.2', fc='#f6e6d8', ec='#9d6f46'))
    for _, row in net_edges.iterrows():
        if row['gene_symbol'] not in gene_pos or row['pathway'] not in path_pos:
            continue
        x1, y1 = gene_pos[row['gene_symbol']]
        x2, y2 = path_pos[row['pathway']]
        ax.plot([x1 + 0.05, x2 - 0.05], [y1, y2], color='#6d6d6d', alpha=0.2 + min(abs(row['NES']) / 4.0, 0.7), linewidth=0.5 + abs(row['NES']) * 0.35)
    save_figure(fig, 'Figure4_gene_pathway_network_optional')


def plot_immune(score_df, expr):
    immune_comp = read_csv(RESULT_DIR / 'tables' / 'immune_primary_comparison.csv')
    top_cells = immune_comp.sort_values(['padj', 'pvalue']).head(6)['cell_type'].tolist()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), sharey=False)
    axes = axes.flatten()
    for ax, cell_type in zip(axes, top_cells):
        sub = score_df[score_df['cell_type'] == cell_type].copy()
        draw_boxplot(ax, sub, 'score', display_cell_label(cell_type), 'Immune ssGSEA NES')
    for ax in axes[len(top_cells):]:
        ax.axis('off')
    save_figure(fig, 'Figure5_immune_score_comparison')

    cor = read_csv(RESULT_DIR / 'tables' / 'gene_immune_correlations.csv')
    if cor.empty:
        return
    pivot = cor.pivot(index='gene_symbol', columns='cell_type', values='rho').fillna(0.0)
    ordered_cols = immune_comp.sort_values(['padj', 'pvalue']).head(8)['cell_type'].tolist()
    pivot = pivot.reindex(columns=[col for col in ordered_cols if col in pivot.columns])
    pivot = pivot.rename(columns={col: display_cell_label(col) for col in pivot.columns})
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.55 * len(pivot))))
    im = draw_heatmap(ax, pivot, 'Gene-immune correlation (Spearman rho)')
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha='right')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    save_figure(fig, 'Figure5_gene_immune_corr_heatmap')

    strongest = cor.sort_values(['padj', 'pvalue', 'rho'], ascending=[True, True, False]).head(4)
    if strongest.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 9.2), constrained_layout=True)
    axes = axes.flatten()
    pivot_score = score_df.pivot(index='sample_id', columns='cell_type', values='score')
    for ax, (_, row) in zip(axes, strongest.iterrows()):
        gene_symbol = row['gene_symbol']
        cell_type = row['cell_type']
        cell_label = display_cell_label(cell_type)
        xs = expr[gene_symbol].loc[pivot_score.index].to_numpy(dtype=float)
        ys = pivot_score[cell_type].to_numpy(dtype=float)
        ax.scatter(xs, ys, s=24, color='#2f6fb0', alpha=0.8)
        coef = np.polyfit(xs, ys, 1)
        xx = np.linspace(xs.min(), xs.max(), 50)
        ax.plot(xx, coef[0] * xx + coef[1], color='#c23b4d', linewidth=1.2)
        ax.set_title(f'{gene_symbol} vs {cell_label}', fontsize=10, pad=8)
        ax.set_xlabel(f'{gene_symbol} log2CPM', fontsize=8, labelpad=4)
        ax.set_ylabel(f'{cell_label} NES', fontsize=8, labelpad=4)
        ax.tick_params(axis='both', labelsize=8)
        ax.grid(alpha=0.18)
    for ax in axes[len(strongest):]:
        ax.axis('off')
    save_figure(fig, 'Supp_gene_immune_scatter')


def plot_qc(pheno):
    qc = read_csv(RESULT_DIR / 'tables' / 'qc_metrics.csv')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    draw_boxplot(axes[0], qc, 'library_size', 'Library size', 'Counts')
    draw_boxplot(axes[1], qc, 'expr_median', 'Expression median', 'Median log2CPM')
    save_figure(fig, 'Supp_QC')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    draw_boxplot(axes[0], qc, 'library_size', 'Library size', 'Counts')
    draw_boxplot(axes[1], qc, 'expr_median', 'Expression median', 'Median log2CPM')
    save_figure(fig, ['library_size', 'expression_distribution'])


def plot_clustering(pheno):
    corr = read_csv(RESULT_DIR / 'tables' / 'sample_correlation_matrix.csv').set_index('sample_id')
    corr = corr.astype(float)
    distance = 1.0 - corr
    distance = (distance + distance.T) / 2.0
    np.fill_diagonal(distance.values, 0.0)
    linkage_mat = linkage(squareform(distance.values, checks=False), method='average')
    order = leaves_list(linkage_mat)
    ordered = corr.iloc[order, order]

    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(ordered.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)
    ax.set_xticks([])
    ax.set_yticks(np.arange(len(ordered.index)))
    ax.set_yticklabels(ordered.index, fontsize=6)
    ax.set_title('Sample correlation heatmap')
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    save_figure(fig, 'sample_correlation_heatmap')
    save_figure(fig, 'Supp_clustering')

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    dendrogram(linkage_mat, labels=ordered.index.tolist(), leaf_rotation=90, leaf_font_size=6, ax=ax)
    ax.set_title('Hierarchical clustering of samples')
    ax.set_ylabel('Distance (1 - correlation)')
    save_figure(fig, 'hierarchical_clustering')


def plot_supplementary_heatmaps(expr, pheno, core_df, deg_df):
    expanded = deg_df[deg_df['significant'] == 1].dropna(subset=['gene_symbol']).head(50)['gene_symbol'].tolist()
    plot_gene_heatmap(expr, pheno, expanded, 'Supp_expanded_deg_heatmap', 'Top significant DE genes')
    plot_gene_heatmap(expr, pheno, core_df['gene_symbol'].tolist(), 'Supp_all_inflammatory_core_heatmap', 'Inflammatory core genes')

    gsea = read_csv(RESULT_DIR / 'tables' / 'gsea_summary_matrix.csv')
    coef_df = read_csv(RESULT_DIR / 'tables' / 'lasso_coefficients.csv')
    genes = coef_df['gene_symbol'].head(4).tolist()
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes = axes.flatten()
    for ax, gene_symbol in zip(axes, genes):
        sub = gsea[gsea['gene_symbol'] == gene_symbol].sort_values(['significant', 'padj', 'NES'], ascending=[False, True, False]).head(6)
        ax.barh(sub['pathway'].str.replace('HALLMARK_', ''), sub['NES'], color=['#c23b4d' if x > 0 else '#2f6fb0' for x in sub['NES']])
        ax.set_title(gene_symbol)
        ax.set_xlabel('NES')
    for ax in axes[len(genes):]:
        ax.axis('off')
    save_figure(fig, 'Supp_single_gene_gsea')


def main():
    ensure_dirs(FIG_DIR, RESULT_DIR / 'logs')
    pheno = read_csv(PROC_DIR / 'pheno_macula_4groups.csv').set_index('sample_id')
    expr = load_symbol_expression()
    score_df = read_csv(RESULT_DIR / 'tables' / 'immune_ssgsea_scores.csv')
    deg_df = read_csv(RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv')
    core_df = read_csv(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv')

    plot_workflow()
    plot_sample_composition(pheno.reset_index())
    plot_pca()
    plot_inflammation()
    plot_volcano()
    plot_gene_heatmap(expr, pheno, core_df['gene_symbol'].tolist(), 'Figure2_deg_inflammatory_heatmap', 'Inflammatory core gene expression')
    plot_overlap(core_df, deg_df)
    plot_lasso_path_and_cv()
    plot_selected_gene_expression(expr, pheno)
    plot_roc()
    plot_gsea()
    plot_immune(score_df, expr)
    plot_qc(pheno)
    plot_clustering(pheno)
    plot_supplementary_heatmaps(expr, pheno, core_df, deg_df)

    log_message('10_make_figures', 'generated matplotlib figures from current result tables')


if __name__ == '__main__':
    main()
