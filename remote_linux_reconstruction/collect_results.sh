#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="${PACKAGE_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
WORK_ROOT="${WORK_ROOT:-${PACKAGE_ROOT}/remote_work}"
ARCHIVE="${WORK_ROOT}/DRnet_GSE276892_remote_results.tar.gz"
HASH_FILE="${ARCHIVE}.sha256"

tar \
    --exclude='*.bam' \
    --exclude='star_index' \
    --exclude='reference' \
    --exclude='tools' \
    -czf "${ARCHIVE}.partial" \
    -C "${PACKAGE_ROOT}" \
    remote_linux_reconstruction/config \
    -C "${WORK_ROOT}" \
    provenance logs \
    all_regions_source_reconstruction/align \
    all_regions_source_reconstruction/lane_counts \
    all_regions_source_reconstruction/results \
    all_regions_source_reconstruction/output_sha256.txt \
    $(if [[ -d "${WORK_ROOT}/primary_assembly_sensitivity/results" ]]; then
        printf '%s ' \
            primary_assembly_sensitivity/align \
            primary_assembly_sensitivity/lane_counts \
            primary_assembly_sensitivity/results \
            primary_assembly_sensitivity/output_sha256.txt
    fi)
mv "${ARCHIVE}.partial" "${ARCHIVE}"
sha256sum "${ARCHIVE}" > "${HASH_FILE}"
echo "Created ${ARCHIVE}"
cat "${HASH_FILE}"
