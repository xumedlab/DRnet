from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


PACKAGE_DIR = Path(__file__).resolve().parents[1]
RESULT_TABLE_DIR = PACKAGE_DIR / "analysis_data" / "results_tables"
SUPP_TABLE_DIR = PACKAGE_DIR / "supplementary_tables"
MAIN_TABLE_DIR = PACKAGE_DIR / "tables"
SUPP_XLSX = PACKAGE_DIR / "Scientific_Reports_supplementary_tables.xlsx"

SEED = 202501
N_PERMUTATIONS = 10000


def main() -> None:
    scores_path = RESULT_TABLE_DIR / "signature_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"Missing required score table: {scores_path}")

    scores = pd.read_csv(scores_path)
    required = {"label", "oof_probability"}
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError(f"Missing required columns in {scores_path}: {sorted(missing)}")

    y = (scores["label"].astype(str) != "healthy control").astype(int).to_numpy()
    score = scores["oof_probability"].astype(float).to_numpy()
    observed_auc = float(roc_auc_score(y, score))

    rng = np.random.default_rng(SEED)
    permuted_aucs = np.empty(N_PERMUTATIONS, dtype=float)
    for i in range(N_PERMUTATIONS):
        permuted_y = rng.permutation(y)
        permuted_aucs[i] = roc_auc_score(permuted_y, score)

    empirical_p = float((np.sum(permuted_aucs >= observed_auc) + 1) / (N_PERMUTATIONS + 1))
    summary = pd.DataFrame(
        [
            {
                "analysis": "fixed_score_label_permutation",
                "score_column": "oof_probability",
                "n_samples": int(len(y)),
                "n_controls": int((y == 0).sum()),
                "n_cases": int((y == 1).sum()),
                "observed_auc": observed_auc,
                "n_permutations": N_PERMUTATIONS,
                "empirical_p_value": empirical_p,
                "null_auc_mean": float(np.mean(permuted_aucs)),
                "null_auc_sd": float(np.std(permuted_aucs, ddof=1)),
                "null_auc_q025": float(np.quantile(permuted_aucs, 0.025)),
                "null_auc_median": float(np.quantile(permuted_aucs, 0.5)),
                "null_auc_q975": float(np.quantile(permuted_aucs, 0.975)),
                "interpretation_boundary": (
                    "Permutation was applied to fixed out-of-fold scores and therefore tests "
                    "label-score separation only; it is not an independent validation and does "
                    "not remove discovery-set feature-selection bias."
                ),
            }
        ]
    )

    SUPP_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    MAIN_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUPP_TABLE_DIR / "Supplementary_Table_S14_label_permutation_check.csv", index=False)
    summary.to_csv(MAIN_TABLE_DIR / "scirep_label_permutation_check.csv", index=False)

    if SUPP_XLSX.exists():
        with pd.ExcelWriter(SUPP_XLSX, mode="a", engine="openpyxl", if_sheet_exists="replace") as writer:
            summary.to_excel(writer, sheet_name="S14_label_permutation", index=False)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
