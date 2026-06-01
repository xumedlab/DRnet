import csv
import runpy

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from pipeline_utils import (
    bh_adjust,
    ensure_dirs,
    load_ensembl_symbol_mapping,
    log_message,
    normalize_ensembl_id,
    spearman,
)

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']

CANDIDATE_COVARIATES = [
    {'column': 'age', 'label': 'age', 'kind': 'continuous'},
    {'column': 'sex', 'label': 'sex', 'kind': 'categorical'},
    {'column': 'post_mortem_interval_min', 'label': 'PMI', 'kind': 'continuous'},
    {'column': 'rin', 'label': 'RIN', 'kind': 'continuous'},
    {'column': 'batch', 'label': 'batch', 'kind': 'categorical'},
]


def load_count_matrix(path):
    df = pd.read_csv(path, sep='\t')
    first_col = df.columns[0]
    if first_col != 'gene':
        df = df.rename(columns={first_col: 'gene'})
    df['gene'] = df['gene'].map(normalize_ensembl_id)
    return df.set_index('gene').T


def prepare_design_metadata(pheno, control_label, case_label):
    subset_meta = pheno[pheno['disease_group'].isin([control_label, case_label])].copy()
    subset_meta['condition'] = pd.Categorical(
        subset_meta['disease_group'],
        categories=[control_label, case_label],
        ordered=True,
    )
    included = []
    design_records = []

    for covariate in CANDIDATE_COVARIATES:
        column = covariate['column']
        label = covariate['label']
        kind = covariate['kind']
        record = {
            'contrast': f'{case_label} vs {control_label}',
            'covariate': label,
            'column': column,
            'kind': kind,
            'status': 'excluded',
            'reason': '',
        }

        if column not in subset_meta.columns:
            record['reason'] = 'not present in metadata'
            design_records.append(record)
            continue

        if kind == 'continuous':
            values = pd.to_numeric(subset_meta[column].replace('', np.nan), errors='coerce')
            missing = int(values.isna().sum())
            scale = float(values.std(ddof=0)) if not values.isna().any() else np.nan
            if missing:
                record['reason'] = f'missing values in {missing} samples'
            elif values.nunique(dropna=True) < 2 or not np.isfinite(scale) or scale == 0:
                record['reason'] = 'no variation in this contrast'
            else:
                subset_meta[column] = (values.astype(float) - float(values.mean())) / scale
                included.append(column)
                record['status'] = 'included'
                record['reason'] = 'complete and variable; z-scored within contrast'
            design_records.append(record)
            continue

        values = subset_meta[column].astype(str).str.strip().replace({'nan': '', 'None': ''})
        missing_mask = values.isin(['', 'NA', 'N/A'])
        missing = int(missing_mask.sum())
        if missing:
            record['reason'] = f'missing values in {missing} samples'
        elif values.nunique(dropna=True) < 2:
            record['reason'] = 'no variation in this contrast'
        elif values.nunique(dropna=True) >= len(values):
            record['reason'] = 'sample-specific or near sample-specific levels'
        else:
            crosstab = pd.crosstab(values, subset_meta['condition'])
            levels_per_condition = (crosstab > 0).sum(axis=1)
            if len(crosstab) <= subset_meta['condition'].nunique() and (levels_per_condition == 1).all():
                record['reason'] = 'fully confounded with condition'
            else:
                subset_meta[column] = pd.Categorical(values)
                included.append(column)
                record['status'] = 'included'
                record['reason'] = 'complete and variable'
        design_records.append(record)

    design_terms = [*included, 'condition']
    design_formula = '~' + ' + '.join(design_terms)
    model_meta = subset_meta[[*included, 'condition']].copy()
    return model_meta, design_formula, included, design_records


def run_deseq_contrast(counts, pheno, control_label, case_label, mapping, out_path):
    subset_meta, design_formula, covariates, design_records = prepare_design_metadata(
        pheno,
        control_label,
        case_label,
    )
    subset_counts = counts.loc[subset_meta.index].astype(int)

    dds = DeseqDataSet(
        counts=subset_counts,
        metadata=subset_meta,
        design=design_formula,
        n_cpus=1,
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=['condition', case_label, control_label], quiet=True, n_cpus=1)
    stat.summary()
    results = stat.results_df.reset_index().rename(
        columns={
            'index': 'gene',
            'log2FoldChange': 'log2FC',
        }
    )

    results['gene'] = results['gene'].map(normalize_ensembl_id)
    results['gene_symbol'] = results['gene'].map(mapping)
    for column in ['baseMean', 'log2FC', 'lfcSE', 'stat', 'pvalue', 'padj']:
        results[column] = pd.to_numeric(results[column], errors='coerce')
    results['pvalue'] = results['pvalue'].fillna(1.0)
    results['padj'] = results['padj'].fillna(1.0)
    results['significant'] = (
        (results['padj'] < cfg['PADJ_THRESHOLD']) & (results['log2FC'].abs() >= cfg['LOG2FC_THRESHOLD'])
    ).astype(int)
    results['significance_rule'] = (
        f"deseq2(padj<{cfg['PADJ_THRESHOLD']},|log2FC|>={cfg['LOG2FC_THRESHOLD']})"
    )
    results['deseq_design_formula'] = design_formula
    results['deseq_covariates'] = ';'.join(covariates) if covariates else 'none'
    results = results.sort_values(['padj', 'pvalue', 'log2FC'], ascending=[True, True, False])
    results.to_csv(out_path, index=False)
    log_message(
        '04_differential_expression',
        (
            f'{control_label} vs {case_label}: n_samples={len(subset_meta)} '
            f'design="{design_formula}" significant={int(results["significant"].sum())}'
        ),
    )
    return design_records


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')

    counts = load_count_matrix(PROC_DIR / 'counts_macula_4groups.tsv')
    log2cpm = load_count_matrix(PROC_DIR / 'log2cpm_macula_4groups.tsv')
    pheno = pd.read_csv(PROC_DIR / 'pheno_macula_4groups.csv').set_index('sample_id')
    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)

    design_records = []
    design_records.extend(run_deseq_contrast(
        counts,
        pheno,
        'healthy control',
        'NPDR/PDR + DME',
        mapping,
        RESULT_DIR / 'tables' / 'deg_primary_healthy_vs_npdr_pdr_dme.csv',
    ))
    design_records.extend(run_deseq_contrast(
        counts,
        pheno,
        'healthy control',
        'diabetic',
        mapping,
        RESULT_DIR / 'tables' / 'deg_diabetic_vs_healthy.csv',
    ))
    design_records.extend(run_deseq_contrast(
        counts,
        pheno,
        'healthy control',
        'NPDR',
        mapping,
        RESULT_DIR / 'tables' / 'deg_nppr_vs_healthy.csv',
    ))
    design_records.extend(run_deseq_contrast(
        counts,
        pheno,
        'diabetic',
        'NPDR/PDR + DME',
        mapping,
        RESULT_DIR / 'tables' / 'deg_npdr_pdr_dme_vs_diabetic.csv',
    ))

    pd.DataFrame(design_records).to_csv(
        RESULT_DIR / 'tables' / 'deseq_covariate_design_summary.csv',
        index=False,
    )

    severity = pheno.loc[log2cpm.index, 'severity_code'].astype(int).tolist()
    trend_rows = []
    for gene in log2cpm.columns:
        values = log2cpm[gene].tolist()
        rho, pvalue = spearman(severity, values)
        trend_rows.append(
            {
                'gene': normalize_ensembl_id(gene),
                'gene_symbol': mapping.get(normalize_ensembl_id(gene), ''),
                'spearman_rho': rho,
                'pvalue': pvalue,
            }
        )

    padj = bh_adjust([row['pvalue'] for row in trend_rows])
    for row, adj in zip(trend_rows, padj):
        row['padj'] = adj
        row['trend_significant'] = int((adj < 0.1) and (row['spearman_rho'] > 0))
        row['trend_rule'] = 'spearman_positive_fdr<0.1'

    with open(RESULT_DIR / 'tables' / 'severity_trend_all_genes.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=trend_rows[0].keys())
        writer.writeheader()
        writer.writerows(sorted(trend_rows, key=lambda x: x['padj']))

    log_message(
        '04_differential_expression',
        f'severity_trend_significant={sum(row["trend_significant"] for row in trend_rows)}',
    )


if __name__ == '__main__':
    main()
