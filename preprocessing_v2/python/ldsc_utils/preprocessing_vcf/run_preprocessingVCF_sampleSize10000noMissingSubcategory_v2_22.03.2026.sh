#!/usr/bin/env bash

set -u

METADATA="/home/mateusz/projects/ifpan-michkor-gr/data/ieu_open_gwas_project/gwas_sampleSize10000_withoutMissingSubcategory/metadata_filtered.tsv"
SCRIPT="/home/mateusz/projects/ifpan-michkor-gr/preprocessing_v2/python/ldsc_utils/ieu_vcf_to_ldsc_and_filter_v2.py"

BASE_OUT="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"
INPUTS_DIR="${BASE_OUT}/inputs"
LOG_DIR="${BASE_OUT}/logs"
SUMMARY_FILE="${BASE_OUT}/preprocessing_summary.tsv"

mkdir -p "${INPUTS_DIR}" "${LOG_DIR}"

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: metadata file not found:"
    echo "${METADATA}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: python script not found:"
    echo "${SCRIPT}"
    exit 1
fi

TOTAL=$(( $(wc -l < "${METADATA}") - 1 ))

if [[ "${TOTAL}" -le 0 ]]; then
    echo "ERROR: metadata file does not contain data rows."
    exit 1
fi

START_TIME=$(date +%s)
DONE=0
OK_COUNT=0
SKIP_COUNT=0
FAIL_COUNT=0

echo -e "id\tstatus\tjob_time_sec\tvcf_file\tout_file\tlog_file" > "${SUMMARY_FILE}"

format_time() {
    local total_sec=$1
    printf '%02d:%02d:%02d' $((total_sec/3600)) $((total_sec%3600/60)) $((total_sec%60))
}

echo "============================================================"
echo "LDSC preprocessing started"
echo "Metadata   : ${METADATA}"
echo "Python     : ${SCRIPT}"
echo "Output dir : ${BASE_OUT}"
echo "Total jobs : ${TOTAL}"
echo "Started at : $(date)"
echo "============================================================"

tail -n +2 "${METADATA}" | while IFS=$'\t' read -r \
    id trait coverage ncase group_name year mr author sex qc_prior_to_upload pmid priority population unit nsnp sample_size build ncontrol covariates subcategory category ontology doi note consortium sd study_design clean_trait folder_name output_dir
do
    DONE=$((DONE + 1))
    JOB_START=$(date +%s)

    vcf_file="${output_dir}/${id}.vcf.gz"
    out_file="${INPUTS_DIR}/${id}.tsv"
    log_file="${LOG_DIR}/${id}.log"

    PERCENT=$(awk "BEGIN {printf \"%.2f\", (${DONE}/${TOTAL})*100}")
    REMAINING_JOBS=$((TOTAL - DONE + 1))

    echo
    echo "------------------------------------------------------------"
    echo "[${DONE}/${TOTAL}] ${PERCENT}% | ${id}"
    echo "VCF : ${vcf_file}"
    echo "OUT : ${out_file}"
    echo "LOG : ${log_file}"
    echo "------------------------------------------------------------"

    if [[ ! -f "${vcf_file}" ]]; then
        JOB_END=$(date +%s)
        JOB_TIME=$((JOB_END - JOB_START))
        SKIP_COUNT=$((SKIP_COUNT + 1))

        echo "[SKIP] Missing VCF file"
        echo -e "${id}\tSKIP_missing_vcf\t${JOB_TIME}\t${vcf_file}\t${out_file}\t${log_file}" >> "${SUMMARY_FILE}"

    else
        if [[ -z "${sample_size}" || "${sample_size}" == "NA" ]]; then
            JOB_END=$(date +%s)
            JOB_TIME=$((JOB_END - JOB_START))
            SKIP_COUNT=$((SKIP_COUNT + 1))

            echo "[SKIP] Missing sample_size"
            echo -e "${id}\tSKIP_missing_sample_size\t${JOB_TIME}\t${vcf_file}\t${out_file}\t${log_file}" >> "${SUMMARY_FILE}"

        else
            stdbuf -oL -eL python3 "${SCRIPT}" \
                --vcf "${vcf_file}" \
                --out "${out_file}" \
                --N "${sample_size}" \
                2>&1 | tee "${log_file}"

            EXIT_CODE=${PIPESTATUS[0]}

            JOB_END=$(date +%s)
            JOB_TIME=$((JOB_END - JOB_START))

            if [[ "${EXIT_CODE}" -eq 0 ]]; then
                OK_COUNT=$((OK_COUNT + 1))
                STATUS="OK"
                echo "[OK] ${id} completed in ${JOB_TIME}s"
            else
                FAIL_COUNT=$((FAIL_COUNT + 1))
                STATUS="FAIL"
                echo "[FAIL] ${id} failed in ${JOB_TIME}s"
            fi

            echo -e "${id}\t${STATUS}\t${JOB_TIME}\t${vcf_file}\t${out_file}\t${log_file}" >> "${SUMMARY_FILE}"
        fi
    fi

    ELAPSED=$(( $(date +%s) - START_TIME ))

    PROCESSED_FOR_AVG=$((OK_COUNT + FAIL_COUNT + SKIP_COUNT))
    if [[ "${PROCESSED_FOR_AVG}" -gt 0 ]]; then
        AVG_TIME=$((ELAPSED / PROCESSED_FOR_AVG))
    else
        AVG_TIME=0
    fi

    REMAINING=$((TOTAL - DONE))
    ETA_SEC=$((AVG_TIME * REMAINING))

    ELAPSED_FMT=$(format_time "${ELAPSED}")
    ETA_FMT=$(format_time "${ETA_SEC}")

    echo "Progress : ${DONE}/${TOTAL} (${PERCENT}%)"
    echo "Counts   : OK=${OK_COUNT} | FAIL=${FAIL_COUNT} | SKIP=${SKIP_COUNT}"
    echo "Last job : ${JOB_TIME}s"
    echo "Average  : ${AVG_TIME}s/job"
    echo "Elapsed  : ${ELAPSED_FMT}"
    echo "ETA      : ${ETA_FMT}"
    echo "Left     : ${REMAINING}"
done

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
TOTAL_TIME_FMT=$(format_time "${TOTAL_TIME}")

echo
echo "============================================================"
echo "LDSC preprocessing finished"
echo "Finished at : $(date)"
echo "Total time  : ${TOTAL_TIME_FMT}"
echo "Summary     : ${SUMMARY_FILE}"
echo "Inputs dir  : ${INPUTS_DIR}"
echo "Logs dir    : ${LOG_DIR}"
echo "============================================================"
