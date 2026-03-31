#!/usr/bin/env bash

set -u
set -o pipefail

############################################
# DEFAULT CONFIG
############################################

BASE_DIR="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"
MUNGED_DIR="${BASE_DIR}/munged"

LDSC="/home/mateusz/projects/ifpan-michkor-gr/tools/ldsc/ldsc.py"
REF="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"
WLD="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"

OUT_DIR="${BASE_DIR}/rg_results_all_vs_all_parallel"
LOG_DIR="${BASE_DIR}/rg_logs_all_vs_all_parallel"
STATUS_DIR="${BASE_DIR}/rg_status_all_vs_all_parallel"
TMP_DIR="${BASE_DIR}/rg_tmp_all_vs_all_parallel"

JOBS=10
SLEEP_SECONDS=2

############################################
# ARGUMENTS
############################################

usage() {
    echo "Usage:"
    echo "  bash run_ldsc_rg_all_vs_all_parallel.sh [-j N]"
    echo
    echo "Options:"
    echo "  -j, --jobs N       Number of parallel jobs (default: 10)"
    echo "  -h, --help         Show this help message"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: missing value for $1"
                usage
                exit 1
            fi
            JOBS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1"
            usage
            exit 1
            ;;
    esac
done

if ! [[ "$JOBS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --jobs must be a positive integer"
    exit 1
fi

if [[ "$JOBS" -lt 1 ]]; then
    echo "ERROR: --jobs must be >= 1"
    exit 1
fi

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

cleanup() {
    echo
    echo "Stopping running background jobs..."
    jobs -rp | xargs -r kill 2>/dev/null
    wait 2>/dev/null
    echo "All running jobs stopped."
    exit 1
}

count_running_jobs() {
    jobs -rp | wc -l
}

wait_for_slot() {
    while true; do
        local running
        running=$(count_running_jobs)

        if [[ "$running" -lt "$JOBS" ]]; then
            break
        fi

        sleep "$SLEEP_SECONDS"
    done
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

run_pair() {
    local file1="$1"
    local file2="$2"

    local trait1
    local trait2
    local out_prefix
    local stdout_log
    local stderr_log
    local status_log
    local start_one
    local end_one
    local runtime
    local exit_code

    trait1=$(basename "$file1" .sumstats.gz)
    trait2=$(basename "$file2" .sumstats.gz)

    out_prefix="${OUT_DIR}/${trait1}_vs_${trait2}"
    stdout_log="${LOG_DIR}/${trait1}_vs_${trait2}.out.log"
    stderr_log="${LOG_DIR}/${trait1}_vs_${trait2}.err.log"
    status_log="${STATUS_DIR}/${trait1}_vs_${trait2}.status.txt"

    if [[ -f "${out_prefix}.log" ]]; then
        {
            echo "status=skipped"
            echo "trait1=${trait1}"
            echo "trait2=${trait2}"
            echo "output_prefix=${out_prefix}"
            echo "stdout_log=${stdout_log}"
            echo "stderr_log=${stderr_log}"
        } > "$status_log"

        echo "SKIP ${trait1} vs ${trait2}"
        return 0
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

    if [[ "$exit_code" -eq 0 ]]; then
        {
            echo "status=ok"
            echo "trait1=${trait1}"
            echo "trait2=${trait2}"
            echo "runtime_sec=${runtime}"
            echo "exit_code=${exit_code}"
            echo "output_prefix=${out_prefix}"
            echo "stdout_log=${stdout_log}"
            echo "stderr_log=${stderr_log}"
        } > "$status_log"

        echo "OK   ${trait1} vs ${trait2} ($(format_seconds "$runtime"))"
    else
        {
            echo "status=failed"
            echo "trait1=${trait1}"
            echo "trait2=${trait2}"
            echo "runtime_sec=${runtime}"
            echo "exit_code=${exit_code}"
            echo "output_prefix=${out_prefix}"
            echo "stdout_log=${stdout_log}"
            echo "stderr_log=${stderr_log}"
        } > "$status_log"

        echo "FAIL ${trait1} vs ${trait2} ($(format_seconds "$runtime"))"
    fi
}

print_header() {
    echo "========================================"
    echo "LDSC all-vs-all parallel run"
    echo "Base dir:        $BASE_DIR"
    echo "Munged dir:      $MUNGED_DIR"
    echo "Results dir:     $OUT_DIR"
    echo "Logs dir:        $LOG_DIR"
    echo "Status dir:      $STATUS_DIR"
    echo "Tmp dir:         $TMP_DIR"
    echo "Traits:          $N"
    echo "Pairs:           $TOTAL_PAIRS"
    echo "Parallel jobs:   $JOBS"
    echo "Threads per job: 1"
    echo "Start time:      $(date)"
    echo "========================================"
    echo
}

print_summary() {
    local ok_count
    local fail_count
    local skip_count
    local done_count
    local total_time

    ok_count=$(find "$STATUS_DIR" -type f -name "*.status.txt" -exec grep -l "^status=ok$" {} + 2>/dev/null | wc -l)
    fail_count=$(find "$STATUS_DIR" -type f -name "*.status.txt" -exec grep -l "^status=failed$" {} + 2>/dev/null | wc -l)
    skip_count=$(find "$STATUS_DIR" -type f -name "*.status.txt" -exec grep -l "^status=skipped$" {} + 2>/dev/null | wc -l)
    done_count=$((ok_count + fail_count + skip_count))
    total_time=$(( $(date +%s) - START_ALL ))

    echo
    echo "========================================"
    echo "RUN FINISHED"
    echo "Traits:          $N"
    echo "Pairs:           $TOTAL_PAIRS"
    echo "Done:            $done_count"
    echo "OK:              $ok_count"
    echo "Failed:          $fail_count"
    echo "Skipped:         $skip_count"
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

trap cleanup INT TERM

START_ALL=$(date +%s)

check_requirements
prepare_dirs
prepare_file_list
prepare_pairs
print_header

while IFS=$'\t' read -r file1 file2; do
    wait_for_slot
    run_pair "$file1" "$file2" &
done < "$PAIR_FILE"

wait
print_summary
