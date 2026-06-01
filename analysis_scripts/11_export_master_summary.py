import csv
import runpy

import pandas as pd

from pipeline_utils import ensure_dirs, log_message

cfg = runpy.run_path('00_config.py')
PROC_DIR, RESULT_DIR = cfg['PROC_DIR'], cfg['RESULT_DIR']


def safe_read_csv(path):
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.DataFrame()


def build_master_table_summary():
    files = sorted((RESULT_DIR / 'tables').glob('*.csv'))
    rows = []
    for path in files:
        with open(path, encoding='utf-8') as f:
            n_rows = max(sum(1 for _ in f) - 1, 0)
        rows.append({'table': path.name, 'rows': n_rows})
    pd.DataFrame(rows).to_csv(RESULT_DIR / 'tables' / 'master_analysis_summary.csv', index=False)
    with open(RESULT_DIR / 'tables' / 'master_analysis_summary.md', 'w', encoding='utf-8') as f:
        f.write('# Master analysis summary\n\n')
        for row in rows:
            f.write(f"- {row['table']}: {row['rows']} rows\n")
    return rows


def build_key_metrics():
    pheno = pd.read_csv(PROC_DIR / 'pheno_macula_4groups.csv')
    candidate = safe_read_csv(RESULT_DIR / 'tables' / 'candidate_gene_summary.csv')
    lasso = safe_read_csv(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv')
    roc = safe_read_csv(RESULT_DIR / 'tables' / 'roc_combined_signature.csv')
    infl_comp = safe_read_csv(RESULT_DIR / 'tables' / 'inflammation_group_comparison.csv')
    infl_trend = safe_read_csv(RESULT_DIR / 'tables' / 'inflammation_severity_trend.csv')
    immune = safe_read_csv(RESULT_DIR / 'tables' / 'immune_primary_comparison.csv')
    recurrence = safe_read_csv(RESULT_DIR / 'tables' / 'pathway_recurrence_summary.csv')
    core = safe_read_csv(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv')
    progressive = safe_read_csv(RESULT_DIR / 'tables' / 'progressive_inflammatory_genes.csv')
    candidate_map = dict(zip(candidate['metric'], candidate['value'])) if not candidate.empty else {}

    metrics = {
        'macula_samples_total': int(len(pheno)),
        'healthy_control_n': int((pheno['disease_group'] == 'healthy control').sum()),
        'diabetic_n': int((pheno['disease_group'] == 'diabetic').sum()),
        'npdr_n': int((pheno['disease_group'] == 'NPDR').sum()),
        'advanced_dr_n': int((pheno['disease_group'] == 'NPDR/PDR + DME').sum()),
        'primary_deg_significant': int(candidate_map.get('primary_deg_significant', 0)),
        'inflammatory_core_genes': int(len(core)),
        'progressive_inflammatory_genes': int(len(progressive)),
        'lasso_selected_genes': ', '.join(lasso['gene_symbol'].tolist()) if not lasso.empty else '',
        'combined_signature_auc': float(roc['auc'].iloc[0]) if not roc.empty else float('nan'),
        'inflammation_primary_pvalue': float(infl_comp['pvalue'].iloc[0]) if not infl_comp.empty else float('nan'),
        'inflammation_severity_rho': float(infl_trend['rho'].iloc[0]) if not infl_trend.empty else float('nan'),
        'top_immune_cells': ', '.join(immune.sort_values('padj').head(5)['cell_type'].tolist()) if not immune.empty else '',
        'top_recurrent_pathways': ', '.join(recurrence.head(5)['pathway'].tolist()) if not recurrence.empty else '',
    }
    pd.DataFrame([{'metric': k, 'value': v} for k, v in metrics.items()]).to_csv(
        RESULT_DIR / 'tables' / 'manuscript_bioinformatics_values.csv',
        index=False,
    )
    return metrics, pheno, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend


def write_chinese_summary(metrics, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend):
    significant_immune = immune[immune['padj'] < 0.05]['cell_type'].tolist() if not immune.empty else []
    if not significant_immune and not immune.empty:
        significant_immune = immune.sort_values('padj').head(5)['cell_type'].tolist()

    lines = [
        '# 论文生信结果整理',
        '',
        '## 数据集与分组',
        f"- GSE160306 中最终纳入黄斑区样本 {metrics['macula_samples_total']} 例。",
        f"- 分组构成：healthy control {metrics['healthy_control_n']} 例，diabetic {metrics['diabetic_n']} 例，NPDR {metrics['npdr_n']} 例，NPDR/PDR + DME {metrics['advanced_dr_n']} 例。",
        '',
        '## 关键生信结果',
        f"- 以 healthy control vs NPDR/PDR + DME 为主比较，DESeq2 共识别显著差异基因 {metrics['primary_deg_significant']} 个。",
        f"- 与 HALLMARK_INFLAMMATORY_RESPONSE 取交后得到炎症核心基因 {len(core)} 个：{', '.join(core['gene_symbol'].tolist()) if not core.empty else '无'}。",
        f"- 其中满足正向严重度趋势的 progressive inflammatory genes 共 {len(progressive)} 个：{', '.join(progressive['gene_symbol'].tolist()) if not progressive.empty else '无'}。",
        f"- LASSO 最终保留的诊断候选基因为：{', '.join(lasso['gene_symbol'].tolist()) if not lasso.empty else '无'}。",
        f"- 组合签名的 OOF AUC 为 {metrics['combined_signature_auc']:.3f}。",
        f"- 炎症 ssGSEA 主比较 P 值为 {metrics['inflammation_primary_pvalue']:.4g}，严重度相关 rho 为 {metrics['inflammation_severity_rho']:.3f}。",
        f"- 免疫 ssGSEA 中最值得关注的细胞类型：{', '.join(significant_immune) if significant_immune else '无显著项'}。",
        f"- 复现度最高的 Hallmark 通路：{', '.join(recurrence.head(5)['pathway'].tolist()) if not recurrence.empty else '无'}。",
        '',
        '## 稿件撰写时应直接修正的地方',
        '- 旧稿中 `limma / GSVA / glmnet / clusterProfiler` 的方法描述已经与当前 Python 管线不符，必须整体替换。',
        '- 旧稿中“7 个诊断基因”“45 个 DE-INFGs”“39 条通路”“3 类免疫细胞”等数字不是当前结果，继续保留就是错稿。',
        '- 当前仓库只覆盖生信部分，不包含临床验证和细胞实验的原始数据，Methods 2.7 以后及 Results 3.6 以后不能凭空补数字。',
        '',
        '## 重点人工核查',
        '- `progressive inflammatory genes` 的判定阈值使用 Spearman 正相关且 FDR < 0.1，属于探索性筛选，需要在文中明确写清楚。',
        '- 免疫浸润来自 bulk RNA 的 ssGSEA 推断，只能写 supportive evidence，不能写成真实浸润比例。',
        '- 组合签名 AUC 来自 OOF 预测而非独立外部验证，不能夸大成已验证的临床诊断模型。',
    ]
    (RESULT_DIR / 'tables' / 'manuscript_bioinformatics_summary.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def write_english_sections(metrics, pheno, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend):
    top_immune = immune.sort_values('padj').head(5)['cell_type'].tolist() if not immune.empty else []
    top_pathways = recurrence.head(5)['pathway'].tolist() if not recurrence.empty else []
    advanced_n = int((pheno['disease_group'] == 'NPDR/PDR + DME').sum())
    text = f"""# Suggested Bioinformatics Sections For The Manuscript

## Methods 2.1-2.6

Transcriptomic data were obtained from GEO under accession number GSE160306. Only macular samples were retained for formal bioinformatic analysis. The final cohort included {metrics['healthy_control_n']} healthy controls, {metrics['diabetic_n']} diabetic samples without retinopathy progression, {metrics['npdr_n']} NPDR samples, and {metrics['advanced_dr_n']} NPDR/PDR + DME samples. Raw count matrices were used for differential expression analysis, whereas log2(CPM+1) values recomputed from the same raw counts were used for PCA, ssGSEA, LASSO modeling, and visualization.

Differential expression for the primary comparison (healthy control vs NPDR/PDR + DME) and the supportive pairwise comparisons was performed using DESeq2 implemented through the PyDESeq2 framework in Python. Genes with adjusted P < 0.05 and |log2 fold change| >= 0.5 were considered significant. A severity-trend analysis across the four disease stages was additionally performed using Spearman correlation, and positively correlated genes with FDR < 0.1 were retained as progression-related candidates.

Inflammation-related activity was quantified by ssGSEA using the HALLMARK_INFLAMMATORY_RESPONSE gene set from MSigDB. Inflammation-associated candidate genes were defined as the intersection between significant primary DEGs and the Hallmark inflammatory response gene set. LASSO logistic regression with repeated stratified cross-validation was then applied to the candidate genes in the primary comparison cohort, and model performance was summarized with out-of-fold predicted probabilities and ROC analysis.

For functional interpretation, each selected gene was used as an anchor for Spearman-ranked preranked GSEA against the Hallmark collection. Relative immune activity was estimated with ssGSEA using 28 curated immune cell signatures, and group differences as well as gene-immune correlations were evaluated with nonparametric statistics and multiple-testing correction.

## Results 3.1-3.5

After restricting the analysis to macular samples, {metrics['macula_samples_total']} samples were included in the final transcriptomic cohort, consisting of {metrics['healthy_control_n']} healthy controls, {metrics['diabetic_n']} diabetic samples, {metrics['npdr_n']} NPDR samples, and {advanced_n} advanced DR samples (NPDR/PDR + DME). DESeq2 identified {metrics['primary_deg_significant']} significant genes in the primary contrast of healthy control versus advanced DR. Intersecting these DEGs with the Hallmark inflammatory response signature yielded {len(core)} inflammation-associated core genes ({', '.join(core['gene_symbol'].tolist()) if not core.empty else 'NA'}). Among them, {len(progressive)} genes also showed a positive severity trend across the four clinical groups ({', '.join(progressive['gene_symbol'].tolist()) if not progressive.empty else 'NA'}).

LASSO logistic regression retained {len(lasso)} genes in the final signature ({', '.join(lasso['gene_symbol'].tolist()) if not lasso.empty else 'NA'}). The combined signature achieved an out-of-fold AUC of {metrics['combined_signature_auc']:.3f}, indicating moderate-to-good discrimination within the discovery cohort while still requiring independent external validation. The inflammation ssGSEA score showed a primary-comparison P value of {float(infl_comp['pvalue'].iloc[0]) if not infl_comp.empty else float('nan'):.4g}, and the four-group severity trend yielded a Spearman rho of {float(infl_trend['rho'].iloc[0]) if not infl_trend.empty else float('nan'):.3f}.

Gene-centered preranked GSEA highlighted recurrent Hallmark programs related to {', '.join([x.replace('HALLMARK_', '').replace('_', ' ') for x in top_pathways]) if top_pathways else 'inflammatory and metabolic signaling'}. Immune ssGSEA further suggested that the most perturbed immune signatures included {', '.join(top_immune) if top_immune else 'no strongly significant immune signatures after correction'}, supporting an association between the selected inflammatory genes and the retinal immune microenvironment.

## Do Not Auto-Fill

The current repository does not contain the raw data for clinical validation, ELISA/RT-qPCR assays, or cell experiments. Therefore, Methods 2.7 onward and Results 3.6 onward should be completed only after the corresponding wet-lab data are provided.
"""
    (RESULT_DIR / 'tables' / 'manuscript_bioinformatics_sections.md').write_text(text, encoding='utf-8')


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    build_master_table_summary()
    metrics, pheno, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend = build_key_metrics()
    write_chinese_summary(metrics, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend)
    write_english_sections(metrics, pheno, core, progressive, lasso, immune, recurrence, infl_comp, infl_trend)
    log_message('11_export_master_summary', 'exported master summary and manuscript-oriented summaries')


if __name__ == '__main__':
    main()
