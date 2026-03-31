#!/usr/bin/env bash

set -u

LDSC_DIR="/home/mateusz/projects/ifpan-michkor-gr/tools/ldsc"
MUNGE_SCRIPT="${LDSC_DIR}/munge_sumstats.py"
MERGE_ALLELES="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/w_hm3.snplist"

BASE_DIR="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"
INPUTS_DIR="${BASE_DIR}/inputs"
MUNGED_DIR="${BASE_DIR}/munged"
LOG_DIR="${BASE_DIR}/logs_munge"
SUMMARY_FILE="${BASE_DIR}/munge_summary.tsv"

mkdir -p "${MUNGED_DIR}" "${LOG_DIR}"

if [[ ! -f "${MUNGE_SCRIPT}" ]]; then
    echo "ERROR: munge_sumstats.py not found:"
    echo "${MUNGE_SCRIPT}"
    exit 1
fi

if [[ ! -f "${MERGE_ALLELES}" ]]; then
    echo "ERROR: HapMap3 SNP list not found:"
    echo "${MERGE_ALLELES}"
    exit 1
fi

mapfile -t INPUT_FILES < <(find "${INPUTS_DIR}" -maxdepth 1 -type f -name "*.tsv" | sort)

TOTAL=${#INPUT_FILES[@]}

if [[ "${TOTAL}" -eq 0 ]]; then
    echo "ERROR: no input TSV files found in:"
    echo "${INPUTS_DIR}"
    exit 1
fi

START_TIME=$(date +%s)
DONE=0
OK_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

echo -e "id\tstatus\tjob_time_sec\tinput_file\tout_prefix\tlog_file" > "${SUMMARY_FILE}"

format_time() {
    local total_sec=$1
    printf '%02d:%02d:%02d' $((total_sec/3600)) $((total_sec%3600/60)) $((total_sec%60))
}

echo "============================================================"
echo "LDSC munging started"
echo "Munge script : ${MUNGE_SCRIPT}"
echo "Input dir    : ${INPUTS_DIR}"
echo "Output dir   : ${MUNGED_DIR}"
echo "Log dir      : ${LOG_DIR}"
echo "HM3 SNP list : ${MERGE_ALLELES}"
echo "Total jobs   : ${TOTAL}"
echo "Started at   : $(date)"
echo "============================================================"

for input_file in "${INPUT_FILES[@]}"; do
    DONE=$((DONE + 1))
    JOB_START=$(date +%s)

    filename=$(basename "${input_file}")
    id="${filename%.tsv}"

    out_prefix="${MUNGED_DIR}/${id}"
    log_file="${LOG_DIR}/${id}.log"

    PERCENT=$(awk "BEGIN {printf \"%.2f\", (${DONE}/${TOTAL})*100}")

    echo
    echo "------------------------------------------------------------"
    echo "[${DONE}/${TOTAL}] ${PERCENT}% | ${id}"
    echo "INPUT : ${input_file}"
    echo "OUT   : ${out_prefix}"
    echo "LOG   : ${log_file}"
    echo "------------------------------------------------------------"

    if [[ ! -f "${input_file}" ]]; then
        JOB_END=$(date +%s)
        JOB_TIME=$((JOB_END - JOB_START))
        SKIP_COUNT=$((SKIP_COUNT + 1))

        echo "[SKIP] Missing input file"
        echo -e "${id}\tSKIP_missing_input\t${JOB_TIME}\t${input_file}\t${out_prefix}\t${log_file}" >> "${SUMMARY_FILE}"
    else
        stdbuf -oL -eL python "${MUNGE_SCRIPT}" \
            --sumstats "${input_file}" \
            --snp SNP \
            --a1 A1 \
            --a2 A2 \
            --signed-sumstats BETA,0 \
            --p P \
            --N-col N \
            --merge-alleles "${MERGE_ALLELES}" \
            --out "${out_prefix}" \
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

        echo -e "${id}\t${STATUS}\t${JOB_TIME}\t${input_file}\t${out_prefix}\t${log_file}" >> "${SUMMARY_FILE}"
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
echo "LDSC munging finished"
echo "Finished at : $(date)"
echo "Total time  : ${TOTAL_TIME_FMT}"
echo "Summary     : ${SUMMARY_FILE}"
echo "Munged dir  : ${MUNGED_DIR}"
echo "Logs dir    : ${LOG_DIR}"
echo "============================================================"
