#!/usr/bin/env bash

set -u
set -o pipefail

############################################
# CONFIG
############################################

BASE_DIR="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"
MUNGED_DIR="${BASE_DIR}/munged"

LDSC="/home/mateusz/projects/ifpan-michkor-gr/tools/ldsc/ldsc.py"
REF="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"
WLD="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"

OUT_DIR="${BASE_DIR}/rg_results_all_vs_all_single_thread"
LOG_DIR="${BASE_DIR}/rg_logs_all_vs_all_single_thread"
STATUS_DIR="${BASE_DIR}/rg_status_all_vs_all_single_thread"
TMP_DIR="${BASE_DIR}/rg_tmp_all_vs_all_single_thread"

############################################
# FUNCTIONS
############################################

format_seconds() {
    local total_seconds="$1"

    if [[ "$total_seconds" -lt 0 ]]; then
        echo "NA"
        return
    fi

    local h=$(( total_seconds / 3600 ))
    local m=$(( (total_seconds % 3600) / 60 ))
    local s=$(( total_seconds % 60 ))

    if [[ "$h" -gt 0 ]]; then
        echo "${h}h ${m}m ${s}s"
    elif [[ "$m" -gt 0 ]]; then
        echo "${m}m ${s}s"
    else
        echo "${s}s"
    fi
}

check_requirements() {
    if [[ ! -d "$BASE_DIR" ]]; then
        echo "ERROR: BASE_DIR not found:"
        echo "$BASE_DIR"
        exit 1
    fi

    if [[ ! -d "$MUNGED_DIR" ]]; then
        echo "ERROR: MUNGED_DIR not found:"
        echo "$MUNGED_DIR"
        exit 1
    fi

    if [[ ! -f "$LDSC" ]]; then
        echo "ERROR: LDSC script not found:"
        echo "$LDSC"
        exit 1
    fi

    if [[ ! -d "$REF" ]]; then
        echo "ERROR: REF directory not found:"
        echo "$REF"
        exit 1
    fi

    if [[ ! -d "$WLD" ]]; then
        echo "ERROR: WLD directory not found:"
        echo "$WLD"
        exit 1
    fi

    if ! command -v python >/dev/null 2>&1; then
        echo "ERROR: python not found in PATH"
        exit 1
    fi
}

prepare_dirs() {
    mkdir -p "$OUT_DIR"
    mkdir -p "$LOG_DIR"
    mkdir -p "$STATUS_DIR"
    mkdir -p "$TMP_DIR"
}

prepare_file_list() {
    FILES=()
    while IFS= read -r -d '' file; do
        FILES+=("$file")
    done < <(find "$MUNGED_DIR" -maxdepth 1 -type f -name "*.sumstats.gz" -print0 | sort -z)

    N=${#FILES[@]}

    if [[ "$N" -lt 2 ]]; then
        echo "ERROR: fewer than 2 munged files found"
        exit 1
    fi
}

prepare_pairs() {
    PAIR_FILE="${TMP_DIR}/pairs.tsv"
    : > "$PAIR_FILE"

    local i
    local j

    for ((i=0; i<N; i++)); do
        for ((j=i+1; j<N; j++)); do
            printf "%s\t%s\n" "${FILES[$i]}" "${FILES[$j]}" >> "$PAIR_FILE"
        done
    done

    TOTAL_PAIRS=$(wc -l < "$PAIR_FILE")

    if [[ "$TOTAL_PAIRS" -eq 0 ]]; then
        echo "ERROR: no pairs generated"
        exit 1
    fi
}

print_header() {
    echo "========================================"
    echo "LDSC all-vs-all single-thread run"
    echo "Base dir:        $BASE_DIR"
    echo "Munged dir:      $MUNGED_DIR"
    echo "Results dir:     $OUT_DIR"
    echo "Logs dir:        $LOG_DIR"
    echo "Status dir:      $STATUS_DIR"
    echo "Tmp dir:         $TMP_DIR"
    echo "Traits:          $N"
    echo "Pairs:           $TOTAL_PAIRS"
    echo "Parallel jobs:   1"
    echo "Threads per job: 1"
    echo "Start time:      $(date)"
    echo "========================================"
    echo
}

print_summary() {
    local total_time
    total_time=$(( $(date +%s) - START_ALL ))

    echo
    echo "========================================"
    echo "RUN FINISHED"
    echo "Traits:          $N"
    echo "Pairs:           $TOTAL_PAIRS"
    echo "Done:            $DONE_TOTAL"
    echo "OK:              $OK_COUNT"
    echo "Failed:          $FAILED_COUNT"
    echo "Skipped:         $SKIPPED_COUNT"
    echo "Total time:      $(format_seconds "$total_time")"
    echo "End time:        $(date)"
    echo "Results dir:     $OUT_DIR"
    echo "Logs dir:        $LOG_DIR"
    echo "Status dir:      $STATUS_DIR"
    echo "========================================"
}

############################################
# MAIN
############################################

START_ALL=$(date +%s)

check_requirements
prepare_dirs
prepare_file_list
prepare_pairs
print_header

DONE_TOTAL=0
OK_COUNT=0
FAILED_COUNT=0
SKIPPED_COUNT=0
SUM_RUNTIME=0
TIMED_RUNS=0

while IFS=$'\t' read -r file1 file2; do
    trait1=$(basename "$file1" .sumstats.gz)
    trait2=$(basename "$file2" .sumstats.gz)

    out_prefix="${OUT_DIR}/${trait1}_vs_${trait2}"
    stdout_log="${LOG_DIR}/${trait1}_vs_${trait2}.out.log"
    stderr_log="${LOG_DIR}/${trait1}_vs_${trait2}.err.log"
    status_log="${STATUS_DIR}/${trait1}_vs_${trait2}.status.txt"

    if [[ -f "${out_prefix}.log" ]]; then
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        DONE_TOTAL=$((DONE_TOTAL + 1))

        remaining=$((TOTAL_PAIRS - DONE_TOTAL))
        elapsed=$(( $(date +%s) - START_ALL ))

        if [[ "$TIMED_RUNS" -gt 0 ]]; then
            avg_runtime=$((SUM_RUNTIME / TIMED_RUNS))
            eta=$((avg_runtime * remaining))
            avg_fmt=$(format_seconds "$avg_runtime")
            eta_fmt=$(format_seconds "$eta")
        else
            avg_fmt="NA"
            eta_fmt="NA"
        fi

        percent=$(awk "BEGIN {printf \"%.2f\", (${DONE_TOTAL}/${TOTAL_PAIRS})*100}")

        {
            echo "status=skipped"
            echo "trait1=${trait1}"
            echo "trait2=${trait2}"
            echo "output_prefix=${out_prefix}"
            echo "stdout_log=${stdout_log}"
            echo "stderr_log=${stderr_log}"
        } > "$status_log"

        echo "[$DONE_TOTAL/$TOTAL_PAIRS] SKIP   ${trait1} vs ${trait2}"
        echo "  progress: ${percent}%"
        echo "  elapsed:  $(format_seconds "$elapsed")"
        echo "  avg time: ${avg_fmt}"
        echo "  ETA:      ${eta_fmt}"
        echo "  failed:   $FAILED_COUNT"
        echo "  skipped:  $SKIPPED_COUNT"
        echo "----------------------------------------"

        continue
    fi

    start_one=$(date +%s)

    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    python "$LDSC" \
        --rg "${file1},${file2}" \
        --ref-ld-chr "$REF" \
        --w-ld-chr "$WLD" \
        --out "$out_prefix" \
        > "$stdout_log" 2> "$stderr_log"

    exit_code=$?

    end_one=$(date +%s)
    runtime=$((end_one - start_one))

    DONE_TOTAL=$((DONE_TOTAL + 1))
    remaining=$((TOTAL_PAIRS - DONE_TOTAL))
    elapsed=$((end_one - START_ALL))

    if [[ "$exit_code" -eq 0 ]]; then
        status="OK"
        OK_COUNT=$((OK_COUNT + 1))
        SUM_RUNTIME=$((SUM_RUNTIME + runtime))
        TIMED_RUNS=$((TIMED_RUNS + 1))
    else
        status="FAIL"
        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi

    if [[ "$TIMED_RUNS" -gt 0 ]]; then
        avg_runtime=$((SUM_RUNTIME / TIMED_RUNS))
        eta=$((avg_runtime * remaining))
        avg_fmt=$(format_seconds "$avg_runtime")
        eta_fmt=$(format_seconds "$eta")
    else
        avg_fmt="NA"
        eta_fmt="NA"
    fi

    percent=$(awk "BEGIN {printf \"%.2f\", (${DONE_TOTAL}/${TOTAL_PAIRS})*100}")

    {
        echo "status=$(echo "$status" | tr '[:upper:]' '[:lower:]')"
        echo "trait1=${trait1}"
        echo "trait2=${trait2}"
        echo "runtime_sec=${runtime}"
        echo "exit_code=${exit_code}"
        echo "output_prefix=${out_prefix}"
        echo "stdout_log=${stdout_log}"
        echo "stderr_log=${stderr_log}"
    } > "$status_log"

    echo "[$DONE_TOTAL/$TOTAL_PAIRS] ${status}    ${trait1} vs ${trait2}"
    echo "  runtime:  $(format_seconds "$runtime")"
    echo "  progress: ${percent}%"
    echo "  elapsed:  $(format_seconds "$elapsed")"
    echo "  avg time: ${avg_fmt}"
    echo "  ETA:      ${eta_fmt}"
    echo "  failed:   $FAILED_COUNT"
    echo "  skipped:  $SKIPPED_COUNT"
    echo "----------------------------------------"

done < "$PAIR_FILE"

print_summary
