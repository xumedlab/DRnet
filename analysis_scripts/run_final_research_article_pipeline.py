#!/usr/bin/env python3
"""Reproduce the final Research Article analyses and run consistency checks."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    package = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, default=package)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=None,
        help="Discovery-input root. Defaults to PACKAGE_ROOT/project_inputs in the release ZIP.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use reduced resampling for a smoke test; outputs are not submission results.",
    )
    return parser.parse_args()


def run(command: list[str], workdir: Path) -> None:
    print("RUN", " ".join(command), flush=True)
    subprocess.run(command, cwd=workdir, check=True)


def require_inputs(project: Path, package: Path) -> None:
    required = [
        project / "data_processed" / "log2cpm_macula_4groups.tsv",
        project / "data_processed" / "ensembl_to_symbol_mapping.csv",
        project / "data_processed" / "manifest_macula_4groups.csv",
        project / "data_raw" / "hallmark_inflammatory_response.gmt",
        project / "GSE102485_expressed_gene_FPKM.txt.gz",
        project / "GSE102485_series_matrix.txt.gz",
        package / "analysis_data" / "independent_validation" / "P2RX4_VALIDATION_PROTOCOL_LOCALLY_FROZEN.md",
        package
        / "analysis_data"
        / "independent_validation"
        / "remote_results"
        / "DRnet_GSE276892_remote_results.tar.gz",
        package
        / "analysis_data"
        / "independent_validation"
        / "remote_results"
        / "DRnet_GSE276892_remote_results.tar.gz.sha256",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required pipeline inputs:\n" + "\n".join(missing))


def main() -> None:
    args = parse_args()
    package = args.package_root.resolve()
    project = (
        args.project_root.resolve()
        if args.project_root is not None
        else (package / "project_inputs").resolve()
    )
    require_inputs(project, package)
    scripts = package / "analysis_scripts"
    python = sys.executable
    discovery_bootstrap = "100" if args.quick else "2000"
    wild_bootstrap = "99" if args.quick else "4999"
    external_bootstrap = "100" if args.quick else "10000"

    run(
        [
            python,
            str(scripts / "33_final_discovery_statistics.py"),
            "--project-root",
            str(project),
            "--package-root",
            str(package),
            "--bootstrap",
            discovery_bootstrap,
            "--wild-bootstrap",
            wild_bootstrap,
        ],
        package,
    )
    run(
        [
            python,
            str(scripts / "38_raw_count_p2rx4_validation.py"),
            "--package-root",
            str(package),
        ],
        package,
    )
    run(
        [
            python,
            str(scripts / "32_independent_p2rx4_validation.py"),
            "--package-root",
            str(package),
            "--bootstrap",
            external_bootstrap,
        ],
        package,
    )
    run(
        [python, str(scripts / "25_voigt_single_cell_localization.py"), "--package-root", str(package)],
        package,
    )
    run([python, str(scripts / "26_updated_study_design.py")], package)

    if args.quick:
        print(
            "QUICK smoke test completed. Re-run without --quick before using the "
            "generated outputs for submission.",
            flush=True,
        )
    else:
        run([python, str(scripts / "34_validate_final_submission.py")], package)


if __name__ == "__main__":
    main()
