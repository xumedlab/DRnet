import csv
import runpy

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from pipeline_utils import (
    ensure_dirs,
    load_ensembl_symbol_mapping,
    log_message,
    matrix_ensembl_to_symbol,
)

cfg = runpy.run_path('00_config.py')
RAW_DIR, PROC_DIR, RESULT_DIR = cfg['RAW_DIR'], cfg['PROC_DIR'], cfg['RESULT_DIR']
SEED = cfg.get('RANDOM_SEED', 202501)


def read_gene_list(path):
    df = pd.read_csv(path)
    if 'gene_symbol' in df.columns:
        return [x for x in df['gene_symbol'].dropna().astype(str) if x]
    if 'gene' in df.columns:
        return [x for x in df['gene'].dropna().astype(str) if x]
    return []


def build_symbol_expression():
    with open(PROC_DIR / 'log2cpm_macula_4groups.tsv', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        matrix_ensembl = {row[0]: [float(x) for x in row[1:]] for row in reader}

    mapping = load_ensembl_symbol_mapping(RAW_DIR, PROC_DIR)
    matrix_symbol, _ = matrix_ensembl_to_symbol(matrix_ensembl, mapping)
    expr = pd.DataFrame(matrix_symbol, index=header[1:]).astype(float)
    expr.index.name = 'sample_id'
    return expr


def bootstrap_auc_ci(y_true, y_score, n_boot=1000):
    rng = np.random.default_rng(SEED)
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    idx = np.arange(len(y_true))
    aucs = []
    for _ in range(n_boot):
        boot_idx = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y_true[boot_idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[boot_idx], y_score[boot_idx]))
    if not aucs:
        return 0.5, 0.5
    aucs = np.sort(np.asarray(aucs))
    return float(np.quantile(aucs, 0.025)), float(np.quantile(aucs, 0.975))


def build_model(c_value):
    return Pipeline(
        [
            ('scaler', StandardScaler()),
            (
                'clf',
                LogisticRegression(
                    penalty='l1',
                    solver='saga',
                    C=float(c_value),
                    max_iter=20000,
                    random_state=SEED,
                ),
            ),
        ]
    )


def score_cs(X, y, cs, cv):
    rows = []
    for c_value in cs:
        fold_scores = []
        for train_idx, valid_idx in cv.split(X, y):
            model = build_model(c_value)
            model.fit(X[train_idx], y[train_idx])
            prob = model.predict_proba(X[valid_idx])[:, 1]
            fold_scores.append(float(roc_auc_score(y[valid_idx], prob)))
        rows.append(
            {
                'C': float(c_value),
                'lambda': float(1.0 / c_value),
                'mean_auc': float(np.mean(fold_scores)),
                'std_auc': float(np.std(fold_scores, ddof=0)),
                'n_folds': len(fold_scores),
            }
        )
    return sorted(rows, key=lambda x: (-x['mean_auc'], x['lambda']))


def main():
    ensure_dirs(RESULT_DIR / 'tables', RESULT_DIR / 'logs')
    progressive = read_gene_list(RESULT_DIR / 'tables' / 'progressive_inflammatory_genes.csv')
    core = read_gene_list(RESULT_DIR / 'tables' / 'inflammatory_core_genes.csv')
    feature_pool = sorted(set(progressive if len(progressive) >= 6 else core))
    if len(feature_pool) < 3:
        raise ValueError('Not enough candidate genes for LASSO signature construction.')

    expr = build_symbol_expression()
    pheno = pd.read_csv(PROC_DIR / 'pheno_macula_4groups.csv').set_index('sample_id')
    primary_ids = pheno[pheno['disease_group'].isin([cfg['PRIMARY_CTRL'], cfg['PRIMARY_CASE']])].index.tolist()
    available = [gene for gene in feature_pool if gene in expr.columns]
    if len(available) < 3:
        raise ValueError('Too few candidate genes are available in the expression matrix after symbol mapping.')

    X = expr.loc[primary_ids, available].to_numpy(dtype=float)
    y = np.asarray(
        [1 if pheno.loc[sample, 'disease_group'] == cfg['PRIMARY_CASE'] else 0 for sample in primary_ids],
        dtype=int,
    )

    cs = np.logspace(-3, 2, 40)
    inner_cv = RepeatedStratifiedKFold(n_splits=4, n_repeats=10, random_state=SEED)
    outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)

    cv_curve = score_cs(X, y, cs, inner_cv)
    with open(RESULT_DIR / 'tables' / 'lasso_cv_curve.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=cv_curve[0].keys())
        writer.writeheader()
        writer.writerows(cv_curve)

    best_c = cv_curve[0]['C']

    path_rows = []
    for c_value in cs:
        model = build_model(c_value)
        model.fit(X, y)
        coef = model.named_steps['clf'].coef_[0]
        for gene_symbol, coefficient in zip(available, coef):
            path_rows.append(
                {
                    'C': float(c_value),
                    'lambda': float(1.0 / c_value),
                    'gene_symbol': gene_symbol,
                    'coefficient': float(coefficient),
                }
            )
    with open(RESULT_DIR / 'tables' / 'lasso_path_coefficients.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=path_rows[0].keys())
        writer.writeheader()
        writer.writerows(path_rows)

    oof_prob = np.zeros(len(y), dtype=float)
    chosen_cs = []
    for train_idx, test_idx in outer_cv.split(X, y):
        inner_rows = score_cs(X[train_idx], y[train_idx], cs, inner_cv)
        fold_c = inner_rows[0]['C']
        chosen_cs.append(fold_c)
        model = build_model(fold_c)
        model.fit(X[train_idx], y[train_idx])
        oof_prob[test_idx] = model.predict_proba(X[test_idx])[:, 1]

    final_c = float(np.median(chosen_cs)) if chosen_cs else float(best_c)
    final_model = build_model(final_c)
    final_model.fit(X, y)
    full_prob = final_model.predict_proba(X)[:, 1]
    decision = final_model.decision_function(X)
    coef = final_model.named_steps['clf'].coef_[0]

    selected = [(gene, float(weight)) for gene, weight in zip(available, coef) if abs(weight) > 1e-8]
    if not selected:
        selected = sorted(
            [(gene, float(weight)) for gene, weight in zip(available, coef)],
            key=lambda item: -abs(item[1]),
        )[: min(6, len(available))]

    selected = sorted(selected, key=lambda item: -abs(item[1]))
    with open(RESULT_DIR / 'tables' / 'lasso_selected_genes.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['gene_symbol'])
        writer.writeheader()
        for gene_symbol, _ in selected:
            writer.writerow({'gene_symbol': gene_symbol})

    with open(RESULT_DIR / 'tables' / 'lasso_coefficients.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['gene_symbol', 'coefficient'])
        writer.writeheader()
        for gene_symbol, coefficient in selected:
            writer.writerow({'gene_symbol': gene_symbol, 'coefficient': coefficient})

    scores = []
    for idx, sample_id in enumerate(primary_ids):
        scores.append(
            {
                'sample_id': sample_id,
                'label': pheno.loc[sample_id, 'disease_group'],
                'signature_score': float(decision[idx]),
                'oof_probability': float(oof_prob[idx]),
                'full_probability': float(full_prob[idx]),
            }
        )
    with open(RESULT_DIR / 'tables' / 'signature_scores.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=scores[0].keys())
        writer.writeheader()
        writer.writerows(scores)

    scaler = final_model.named_steps['scaler']
    X_scaled = scaler.transform(X)
    roc_rows = []
    roc_curve_rows = []
    for gene_symbol, _ in selected:
        col_idx = available.index(gene_symbol)
        gene_score = X_scaled[:, col_idx]
        auc = float(roc_auc_score(y, gene_score))
        lo, hi = bootstrap_auc_ci(y, gene_score)
        roc_rows.append({'item': gene_symbol, 'auc': auc, 'ci95_low': lo, 'ci95_high': hi})
        fpr, tpr, thresholds = roc_curve(y, gene_score)
        for fp, tp, threshold in zip(fpr, tpr, thresholds):
            roc_curve_rows.append(
                {
                    'item': gene_symbol,
                    'fpr': float(fp),
                    'tpr': float(tp),
                    'threshold': float(threshold),
                }
            )

    with open(RESULT_DIR / 'tables' / 'roc_individual_genes.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=roc_rows[0].keys())
        writer.writeheader()
        writer.writerows(sorted(roc_rows, key=lambda row: -row['auc']))

    with open(RESULT_DIR / 'tables' / 'roc_curve_individual_genes.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=roc_curve_rows[0].keys())
        writer.writeheader()
        writer.writerows(roc_curve_rows)

    combined_auc = float(roc_auc_score(y, oof_prob))
    combined_lo, combined_hi = bootstrap_auc_ci(y, oof_prob)
    with open(RESULT_DIR / 'tables' / 'roc_combined_signature.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['item', 'auc', 'ci95_low', 'ci95_high', 'lambda'])
        writer.writeheader()
        writer.writerow(
            {
                'item': 'combined_signature_oof',
                'auc': combined_auc,
                'ci95_low': combined_lo,
                'ci95_high': combined_hi,
                'lambda': float(1.0 / final_c),
            }
        )

    fpr, tpr, thresholds = roc_curve(y, oof_prob)
    with open(RESULT_DIR / 'tables' / 'roc_curve_combined_signature.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['item', 'fpr', 'tpr', 'threshold'])
        writer.writeheader()
        for fp, tp, threshold in zip(fpr, tpr, thresholds):
            writer.writerow(
                {
                    'item': 'combined_signature_oof',
                    'fpr': float(fp),
                    'tpr': float(tp),
                    'threshold': float(threshold),
                }
            )

    log_message(
        '07_lasso_signature',
        f'feature_pool={len(available)} selected={len(selected)} auc={combined_auc:.4f} lambda={1.0 / final_c:.6f}',
    )


if __name__ == '__main__':
    main()
