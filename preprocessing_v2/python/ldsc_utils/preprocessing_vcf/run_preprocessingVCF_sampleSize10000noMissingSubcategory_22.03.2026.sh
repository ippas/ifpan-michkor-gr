#!/usr/bin/env bash

set -u

METADATA="/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/gwas_sampleSize10000_withoutMissingSubcategory/metadata_filtered.tsv"
SCRIPT="/home/mateusz/projects/ifpan-michkor-gr/preprocessing_v2/python/ldsc_utils/ieu_vcf_to_ldsc_and_filter_v2.py"

BASE_OUT="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"
INPUTS_DIR="${BASE_OUT}/inputs"
LOG_DIR="${BASE_OUT}/logs"

mkdir -p "${INPUTS_DIR}" "${LOG_DIR}"

TOTAL=$(($(wc -l < "$METADATA") - 1))
DONE=0
START_TIME=$(date +%s)

echo "Total jobs: $TOTAL"
echo "----------------------------------------"

tail -n +2 "${METADATA}" | while IFS=$'\t' read -r \
    id trait coverage ncase group_name year mr author sex qc_prior_to_upload pmid priority population unit nsnp sample_size build ncontrol covariates subcategory category ontology doi note consortium sd study_design clean_trait folder_name output_dir
do
    JOB_START=$(date +%s)

    vcf_file="${output_dir}/${id}.vcf.gz"
    out_file="${INPUTS_DIR}/${id}.tsv"
    log_file="${LOG_DIR}/${id}.log"

    if [[ ! -f "${vcf_file}" ]]; then
        echo "[SKIP] ${id} (no VCF)"
        continue
    fi

    if [[ -z "${sample_size}" || "${sample_size}" == "NA" ]]; then
        echo "[SKIP] ${id} (no sample_size)"
        continue
    fi

    python3 "${SCRIPT}" \
        --vcf "${vcf_file}" \
        --out "${out_file}" \
        --N "${sample_size}" \
        > "${log_file}" 2>&1

    JOB_END=$(date +%s)
    JOB_TIME=$((JOB_END - JOB_START))

    DONE=$((DONE + 1))

    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TIME))

    AVG_TIME=$((ELAPSED / DONE))
    REMAINING=$((TOTAL - DONE))
    ETA_SEC=$((AVG_TIME * REMAINING))

    PERCENT=$(awk "BEGIN {printf \"%.2f\", (${DONE}/${TOTAL})*100}")

    # format czasu
    ETA_FMT=$(printf '%02d:%02d:%02d' $((ETA_SEC/3600)) $((ETA_SEC%3600/60)) $((ETA_SEC%60)))
    ELAPSED_FMT=$(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60)))

    echo "[${DONE}/${TOTAL}] (${PERCENT}%) | job=${JOB_TIME}s | avg=${AVG_TIME}s | elapsed=${ELAPSED_FMT} | ETA=${ETA_FMT} | ${id}"

done
