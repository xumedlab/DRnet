#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="${PACKAGE_ROOT:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
WORK_ROOT="${WORK_ROOT:-${PACKAGE_ROOT}/remote_work}"
THREADS="${THREADS:-16}"
KEEP_BAM="${KEEP_BAM:-0}"
RUN_PRIMARY_ASSEMBLY_SENSITIVITY="${RUN_PRIMARY_ASSEMBLY_SENSITIVITY:-1}"
SJDB_OVERHANG="${SJDB_OVERHANG:-100}"
CLEAN_GENERATED_REFERENCE_AND_INDEX="${CLEAN_GENERATED_REFERENCE_AND_INDEX:-0}"

mkdir -p "${WORK_ROOT}/logs" "${WORK_ROOT}/provenance"
RUN_LOG="${WORK_ROOT}/logs/reconstruction_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "Start: $(date --iso-8601=seconds)"
echo "PACKAGE_ROOT=${PACKAGE_ROOT}"
echo "WORK_ROOT=${WORK_ROOT}"
echo "THREADS=${THREADS}"
echo "KEEP_BAM=${KEEP_BAM}"
echo "RUN_PRIMARY_ASSEMBLY_SENSITIVITY=${RUN_PRIMARY_ASSEMBLY_SENSITIVITY}"
echo "SJDB_OVERHANG=${SJDB_OVERHANG}"
echo "CLEAN_GENERATED_REFERENCE_AND_INDEX=${CLEAN_GENERATED_REFERENCE_AND_INDEX}"

bash "${SCRIPT_DIR}/preflight.sh"

STAR_BIN="${WORK_ROOT}/tools/STAR-2.7.8a/bin/Linux_x86_64_static/STAR"
FEATURECOUNTS_BIN="${WORK_ROOT}/tools/subread-2.0.1-Linux-x86_64/bin/featureCounts"
SAMPLE_SHEET="${SCRIPT_DIR}/config/sample_sheet.tsv"

run_workflow() {
    local workflow_label="$1"
    local fasta_gz="$2"
    local gtf_gz="$3"
    local workflow_root="${WORK_ROOT}/${workflow_label}"
    local reference_root="${workflow_root}/reference"
    local index_root="${workflow_root}/star_index"
    local align_root="${workflow_root}/align"
    local counts_root="${workflow_root}/lane_counts"
    local results_root="${workflow_root}/results"
    mkdir -p "${reference_root}" "${align_root}" "${counts_root}" "${results_root}"

    local fasta="${reference_root}/genome.fa"
    local gtf="${reference_root}/annotation.gtf"
    if [[ ! -s "${fasta}" ]]; then
        echo "[$(date --iso-8601=seconds)] Decompressing ${workflow_label} FASTA"
        gzip -dc "${fasta_gz}" > "${fasta}.partial"
        mv "${fasta}.partial" "${fasta}"
    fi
    if [[ ! -s "${gtf}" ]]; then
        echo "[$(date --iso-8601=seconds)] Decompressing ${workflow_label} GTF"
        gzip -dc "${gtf_gz}" > "${gtf}.partial"
        mv "${gtf}.partial" "${gtf}"
    fi

    if [[ ! -s "${index_root}/Genome" ]]; then
        echo "[$(date --iso-8601=seconds)] Building STAR index: ${workflow_label}"
        local index_build="${workflow_root}/star_index_building"
        local index_stamp
        index_stamp="$(date +%Y%m%d_%H%M%S)"
        if [[ -d "${index_build}" ]]; then
            mv "${index_build}" "${index_build}.failed_${index_stamp}"
        fi
        mkdir -p "${index_build}"
        "${STAR_BIN}" \
            --runMode genomeGenerate \
            --runThreadN "${THREADS}" \
            --genomeDir "${index_build}" \
            --genomeFastaFiles "${fasta}" \
            --sjdbGTFfile "${gtf}" \
            --sjdbOverhang "${SJDB_OVERHANG}"
        if [[ ! -s "${index_build}/Genome" ]]; then
            echo "ERROR: STAR index build did not create Genome for ${workflow_label}" >&2
            exit 1
        fi
        if [[ -d "${index_root}" ]]; then
            mv "${index_root}" "${index_root}.incomplete_${index_stamp}"
        fi
        mv "${index_build}" "${index_root}"
    else
        echo "[$(date --iso-8601=seconds)] Reusing complete STAR index: ${workflow_label}"
    fi

    while IFS=$'\t' read -r run_accession sample_id diagnosis source_dataset is_reused_control lane_index lanes_per_sample fastq_relative_path expected_bytes expected_md5 expected_sha256 official_fastq_url sample_accession experiment_accession secondary_sample_accession library_layout instrument_model; do
        if [[ "${run_accession}" == "run_accession" ]]; then
            continue
        fi
        if [[ "${library_layout}" != "SINGLE" ]]; then
            echo "ERROR: unsupported library layout for ${run_accession}: ${library_layout}" >&2
            exit 1
        fi
        local fastq="${PACKAGE_ROOT}/${fastq_relative_path}"
        local run_root="${align_root}/${run_accession}"
        local counts_file="${counts_root}/${run_accession}.featureCounts.txt"
        local counts_summary="${counts_file}.summary"
        mkdir -p "${run_root}"
        if [[ -s "${counts_file}" && -s "${counts_summary}" && -s "${run_root}/Log.final.out" ]]; then
            echo "[$(date --iso-8601=seconds)] Resume: ${workflow_label}/${run_accession} already counted"
            continue
        fi
        echo "[$(date --iso-8601=seconds)] Aligning ${workflow_label}/${run_accession} (${sample_id}, lane ${lane_index}/${lanes_per_sample})"
        "${STAR_BIN}" \
            --runThreadN "${THREADS}" \
            --genomeDir "${index_root}" \
            --readFilesIn "${fastq}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${run_root}/" \
            --outSAMtype BAM Unsorted \
            --quantMode GeneCounts
        local bam="${run_root}/Aligned.out.bam"
        if [[ ! -s "${bam}" || ! -s "${run_root}/Log.final.out" ]]; then
            echo "ERROR: STAR did not create a complete BAM/log for ${run_accession}" >&2
            exit 1
        fi
        echo "[$(date --iso-8601=seconds)] Counting ${workflow_label}/${run_accession}"
        "${FEATURECOUNTS_BIN}" \
            -T "${THREADS}" \
            -a "${gtf}" \
            -o "${counts_file}.partial" \
            -t exon \
            -g gene_id \
            -s 0 \
            "${bam}"
        mv "${counts_file}.partial" "${counts_file}"
        mv "${counts_file}.partial.summary" "${counts_summary}"
        if [[ ! -s "${counts_file}" || ! -s "${counts_summary}" ]]; then
            echo "ERROR: featureCounts did not create complete output for ${run_accession}" >&2
            exit 1
        fi
        if [[ "${KEEP_BAM}" == "0" ]]; then
            rm -f "${bam}"
        fi
    done < "${SAMPLE_SHEET}"

    python3 "${SCRIPT_DIR}/scripts/aggregate_featurecounts.py" \
        --sample-sheet "${SAMPLE_SHEET}" \
        --workflow-root "${workflow_root}" \
        --gtf "${gtf}" \
        --output-dir "${results_root}"

    (
        cd "${workflow_root}"
        find results lane_counts align -type f ! -name '*.bam' -print0 \
            | sort -z \
            | xargs -0 sha256sum
    ) > "${workflow_root}/output_sha256.txt"
    if [[ "${CLEAN_GENERATED_REFERENCE_AND_INDEX}" == "1" ]]; then
        case "${index_root}" in
            "${workflow_root}/"*) ;;
            *)
                echo "ERROR: refusing to clean index outside workflow root: ${index_root}" >&2
                exit 1
                ;;
        esac
        case "${reference_root}" in
            "${workflow_root}/"*) ;;
            *)
                echo "ERROR: refusing to clean reference outside workflow root: ${reference_root}" >&2
                exit 1
                ;;
        esac
        echo "[$(date --iso-8601=seconds)] Removing reproducible index/reference to release disk: ${workflow_label}"
        rm -rf -- "${index_root}" "${reference_root}"
    fi
    echo "[$(date --iso-8601=seconds)] Completed workflow: ${workflow_label}"
}

run_workflow \
    "all_regions_source_reconstruction" \
    "${PACKAGE_ROOT}/tools/reference/gencode_v42_all_regions/GRCh38.p13.genome.fa.gz" \
    "${PACKAGE_ROOT}/tools/reference/gencode_v42_all_regions/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz"

if [[ "${RUN_PRIMARY_ASSEMBLY_SENSITIVITY}" == "1" ]]; then
    run_workflow \
        "primary_assembly_sensitivity" \
        "${PACKAGE_ROOT}/tools/reference/gencode_v42_primary_assembly/GRCh38.primary_assembly.genome.fa.gz" \
        "${PACKAGE_ROOT}/tools/reference/gencode_v42_primary_assembly/gencode.v42.primary_assembly.annotation.gtf.gz"
fi

bash "${SCRIPT_DIR}/collect_results.sh"
echo "Finish: $(date --iso-8601=seconds)"
echo "Result bundle: ${WORK_ROOT}/DRnet_GSE276892_remote_results.tar.gz"
