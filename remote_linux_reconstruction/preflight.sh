#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="${PACKAGE_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
WORK_ROOT="${WORK_ROOT:-${PACKAGE_ROOT}/remote_work}"
MIN_FREE_GB="${MIN_FREE_GB:-180}"

mkdir -p "${WORK_ROOT}/provenance" "${WORK_ROOT}/tools"

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "ERROR: this workflow must run on Linux." >&2
    exit 1
fi
if [[ "$(uname -m)" != "x86_64" ]]; then
    echo "ERROR: official bundled binaries require Linux x86_64." >&2
    exit 1
fi
for command_name in bash python3 tar gzip md5sum sha256sum awk sed; do
    command -v "${command_name}" >/dev/null || {
        echo "ERROR: missing command: ${command_name}" >&2
        exit 1
    }
done

available_kb="$(df -Pk "${WORK_ROOT}" | awk 'NR==2 {print $4}')"
available_gb="$((available_kb / 1024 / 1024))"
memory_kb="$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
memory_gb="$((memory_kb / 1024 / 1024))"
echo "Linux architecture: $(uname -m)"
echo "Detected RAM: ${memory_gb} GiB"
echo "Input/package filesystem: $(df -hP "${PACKAGE_ROOT}" | awk 'NR==2 {print $1, $4, $6}')"
echo "Work filesystem: $(df -hP "${WORK_ROOT}" | awk 'NR==2 {print $1, $4, $6}')"
echo "Free disk at WORK_ROOT: ${available_gb} GiB"
if (( memory_gb < 32 )); then
    echo "ERROR: less than 32 GiB RAM; use a server with at least 32 GiB, preferably 64 GiB." >&2
    exit 1
fi
if (( available_gb < MIN_FREE_GB )); then
    echo "ERROR: WORK_ROOT has less than ${MIN_FREE_GB} GiB free disk. Increase capacity or set MIN_FREE_GB only after reviewing peak usage." >&2
    exit 1
fi

python3 "${SCRIPT_DIR}/scripts/verify_remote_inputs.py" \
    --package-root "${PACKAGE_ROOT}" \
    --manifest "${SCRIPT_DIR}/config/input_files.tsv" \
    --report "${WORK_ROOT}/provenance/remote_input_verification.tsv"

STAR_ARCHIVE="${PACKAGE_ROOT}/tools/downloads/STAR-2.7.8a.tar.gz"
SUBREAD_ARCHIVE="${PACKAGE_ROOT}/tools/downloads/subread-2.0.1-Linux-x86_64.tar.gz"
if [[ ! -x "${WORK_ROOT}/tools/STAR-2.7.8a/bin/Linux_x86_64_static/STAR" ]]; then
    tar -xzf "${STAR_ARCHIVE}" -C "${WORK_ROOT}/tools"
    chmod +x "${WORK_ROOT}/tools/STAR-2.7.8a/bin/Linux_x86_64_static/STAR"
fi
if [[ ! -x "${WORK_ROOT}/tools/subread-2.0.1-Linux-x86_64/bin/featureCounts" ]]; then
    tar -xzf "${SUBREAD_ARCHIVE}" -C "${WORK_ROOT}/tools"
    chmod +x "${WORK_ROOT}/tools/subread-2.0.1-Linux-x86_64/bin/featureCounts"
fi

STAR_BIN="${WORK_ROOT}/tools/STAR-2.7.8a/bin/Linux_x86_64_static/STAR"
FEATURECOUNTS_BIN="${WORK_ROOT}/tools/subread-2.0.1-Linux-x86_64/bin/featureCounts"
{
    date --iso-8601=seconds
    uname -a
    echo "STAR: $(${STAR_BIN} --version)"
    "${FEATURECOUNTS_BIN}" -v 2>&1
    python3 --version
    df -h "${PACKAGE_ROOT}" "${WORK_ROOT}"
    free -h
} > "${WORK_ROOT}/provenance/environment_preflight.txt"

echo "Preflight PASS. Report: ${WORK_ROOT}/provenance/remote_input_verification.tsv"
