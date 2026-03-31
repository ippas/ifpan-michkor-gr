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
SLEEP_SECONDS=1

TIMEOUT_SECONDS=300          # 5 minutes
MAX_TIMEOUT_RETRIES=1        # retry once only after timeout
THREADS_PER_JOB=1

############################################
# ARGUMENTS
############################################

usage() {
    echo "Usage:"
    echo "  bash run_ldsc_rg_all_vs_all_parallel_v3.sh [options]"
    echo
    echo "Options:"
    echo "  -j, --jobs N                Number of parallel jobs (default: 10)"
    echo "  -t, --timeout-seconds N     Timeout per pair in seconds (default: 300)"
    echo "  -r, --timeout-retries N     Retries only after timeout (default: 1)"
    echo "  -h, --help                  Show this help message"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)
            [[ $# -ge 2 ]] || { echo "ERROR: missing value for $1"; usage; exit 1; }
            JOBS="$2"
            shift 2
            ;;
        -t|--timeout-seconds)
            [[ $# -ge 2 ]] || { echo "ERROR: missing value for $1"; usage; exit 1; }
            TIMEOUT_SECONDS="$2"
            shift 2
            ;;
        -r|--timeout-retries)
            [[ $# -ge 2 ]] || { echo "ERROR: missing value for $1"; usage; exit 1; }
            MAX_TIMEOUT_RETRIES="$2"
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

[[ "$JOBS" =~ ^[0-9]+$ ]] || { echo "ERROR: --jobs must be a positive integer"; exit 1; }
[[ "$TIMEOUT_SECONDS" =~ ^[0-9]+$ ]] || { echo "ERROR: --timeout-seconds must be a positive integer"; exit 1; }
[[ "$MAX_TIMEOUT_RETRIES" =~ ^[0-9]+$ ]] || { echo "ERROR: --timeout-retries must be >= 0"; exit 1; }

[[ "$JOBS" -ge 1 ]] || { echo "ERROR: --jobs must be >= 1"; exit 1; }
[[ "$TIMEOUT_SECONDS" -ge 1 ]] || { echo "ERROR: --timeout-seconds must be >= 1"; exit 1; }

############################################
# GLOBALS
############################################

PAIR_FILE=""
START_ALL=0
N=0
TOTAL_PAIRS=0

LOCK_FILE=""
DONE_COUNT_FILE=""
OK_COUNT_FILE=""
FAIL_COUNT_FILE=""
SKIP_COUNT_FILE=""
TIMEOUT_COUNT_FILE=""
RETRY_COUNT_FILE=""

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
    [[ -d "$BASE_DIR" ]] || { echo "ERROR: BASE_DIR not found: $BASE_DIR"; exit 1; }
    [[ -d "$MUNGED_DIR" ]] || { echo "ERROR: MUNGED_DIR not found: $MUNGED_DIR"; exit 1; }
    [[ -f "$LDSC" ]] || { echo "ERROR: LDSC script not found: $LDSC"; exit 1; }
    [[ -d "$REF" ]] || { echo "ERROR: REF directory not found: $REF"; exit 1; }
    [[ -d "$WLD" ]] || { echo "ERROR: WLD directory not found: $WLD"; exit 1; }

    command -v python >/dev/null 2>&1 || { echo "ERROR: python not found in PATH"; exit 1; }
    command -v timeout >/dev/null 2>&1 || { echo "ERROR: timeout command not found in PATH"; exit 1; }
    command -v flock >/dev/null 2>&1 || { echo "ERROR: flock command not found in PATH"; exit 1; }
}

prepare_dirs() {
    mkdir -p "$OUT_DIR" "$LOG_DIR" "$STATUS_DIR" "$TMP_DIR"

    LOCK_FILE="${TMP_DIR}/progress.lock"
    DONE_COUNT_FILE="${TMP_DIR}/done.count"
    OK_COUNT_FILE="${TMP_DIR}/ok.count"
    FAIL_COUNT_FILE="${TMP_DIR}/fail.count"
    SKIP_COUNT_FILE="${TMP_DIR}/skip.count"
    TIMEOUT_COUNT_FILE="${TMP_DIR}/timeout.count"
    RETRY_COUNT_FILE="${TMP_DIR}/retry.count"

    echo 0 > "$DONE_COUNT_FILE"
    echo 0 > "$OK_COUNT_FILE"
    echo 0 > "$FAIL_COUNT_FILE"
    echo 0 > "$SKIP_COUNT_FILE"
    echo 0 > "$TIMEOUT_COUNT_FILE"
    echo 0 > "$RETRY_COUNT_FILE"

    : > "$LOCK_FILE"
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

cleanup() {
    echo
    echo "Stopping running background jobs..."
    jobs -rp | xargs -r kill 2>/dev/null || true
    wait 2>/dev/null || true
    echo "All running jobs stopped."
    exit 1
}

counter_add() {
    local file="$1"
    local delta="$2"

    (
        flock -x 200
        local current=0
        [[ -f "$file" ]] && read -r current < "$file"
        echo $(( current + delta )) > "$file"
    ) 200>"$LOCK_FILE"
}

get_counter() {
    local file="$1"
    local value=0
    [[ -f "$file" ]] && read -r value < "$file"
    echo "$value"
}

print_progress_line() {
    local label="$1"
    local trait1="$2"
    local trait2="$3"
    local runtime="$4"
    local attempt="$5"

    (
        flock -x 200

        local done_count ok_count fail_count skip_count timeout_count retry_count left_count
        local percent

        read -r done_count < "$DONE_COUNT_FILE"
        read -r ok_count < "$OK_COUNT_FILE"
        read -r fail_count < "$FAIL_COUNT_FILE"
        read -r skip_count < "$SKIP_COUNT_FILE"
        read -r timeout_count < "$TIMEOUT_COUNT_FILE"
        read -r retry_count < "$RETRY_COUNT_FILE"

        left_count=$(( TOTAL_PAIRS - done_count ))
        percent=$(awk -v d="$done_count" -v t="$TOTAL_PAIRS" 'BEGIN {printf "%.2f", 100*d/t}')

        if [[ "$runtime" -ge 0 ]]; then
            echo "[${done_count}/${TOTAL_PAIRS} | ${percent}% | left:${left_count} | ok:${ok_count} fail:${fail_count} timeout:${timeout_count} skip:${skip_count} retry:${retry_count}] ${label} ${trait1} vs ${trait2} (attempt:${attempt}, $(format_seconds "$runtime"))"
        else
            echo "[${done_count}/${TOTAL_PAIRS} | ${percent}% | left:${left_count} | ok:${ok_count} fail:${fail_count} timeout:${timeout_count} skip:${skip_count} retry:${retry_count}] ${label} ${trait1} vs ${trait2} (attempt:${attempt})"
        fi
    ) 200>"$LOCK_FILE"
}

write_status_file() {
    local status_log="$1"
    local status="$2"
    local trait1="$3"
    local trait2="$4"
    local runtime="$5"
    local exit_code="$6"
    local out_prefix="$7"
    local stdout_log="$8"
    local stderr_log="$9"
    local attempt="${10}"

    {
        echo "status=${status}"
        echo "trait1=${trait1}"
        echo "trait2=${trait2}"
        echo "runtime_sec=${runtime}"
        echo "exit_code=${exit_code}"
        echo "output_prefix=${out_prefix}"
        echo "stdout_log=${stdout_log}"
        echo "stderr_log=${stderr_log}"
        echo "attempt=${attempt}"
        echo "timestamp=$(date '+%F %T')"
    } > "$status_log"
}

run_pair() {
    local file1="$1"
    local file2="$2"

    local trait1 trait2 pair_name
    local out_prefix status_log done_marker
    local start_one end_one runtime
    local attempt max_attempts
    local stdout_log stderr_log
    local exit_code final_status

    trait1=$(basename "$file1" .sumstats.gz)
    trait2=$(basename "$file2" .sumstats.gz)
    pair_name="${trait1}_vs_${trait2}"

    out_prefix="${OUT_DIR}/${pair_name}"
    status_log="${STATUS_DIR}/${pair_name}.status.txt"
    done_marker="${STATUS_DIR}/${pair_name}.done"

    if [[ -f "$done_marker" ]] || { [[ -f "$status_log" ]] && grep -q '^status=ok$' "$status_log"; }; then
        counter_add "$SKIP_COUNT_FILE" 1
        counter_add "$DONE_COUNT_FILE" 1
        print_progress_line "SKIP" "$trait1" "$trait2" 0 0
        return 0
    fi

    max_attempts=$(( MAX_TIMEOUT_RETRIES + 1 ))
    attempt=1
    final_status="failed"

    while [[ "$attempt" -le "$max_attempts" ]]; do
        stdout_log="${LOG_DIR}/${pair_name}.attempt${attempt}.out.log"
        stderr_log="${LOG_DIR}/${pair_name}.attempt${attempt}.err.log"

        start_one=$(date +%s)

        OMP_NUM_THREADS="$THREADS_PER_JOB" \
        OPENBLAS_NUM_THREADS="$THREADS_PER_JOB" \
        MKL_NUM_THREADS="$THREADS_PER_JOB" \
        NUMEXPR_NUM_THREADS="$THREADS_PER_JOB" \
        VECLIB_MAXIMUM_THREADS="$THREADS_PER_JOB" \
        timeout --signal=TERM --kill-after=30s "${TIMEOUT_SECONDS}s" \
        python "$LDSC" \
            --rg "${file1},${file2}" \
            --ref-ld-chr "$REF" \
            --w-ld-chr "$WLD" \
            --out "$out_prefix" \
            > "$stdout_log" 2> "$stderr_log"

        exit_code=$?

        end_one=$(date +%s)
        runtime=$(( end_one - start_one ))

        if [[ "$exit_code" -eq 0 ]]; then
            final_status="ok"
            write_status_file "$status_log" "$final_status" "$trait1" "$trait2" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"
            touch "$done_marker"

            counter_add "$OK_COUNT_FILE" 1
            counter_add "$DONE_COUNT_FILE" 1
            print_progress_line "OK  " "$trait1" "$trait2" "$runtime" "$attempt"
            return 0
        fi

        if [[ "$exit_code" -eq 124 || "$exit_code" -eq 137 ]]; then
            if [[ "$attempt" -lt "$max_attempts" ]]; then
                write_status_file "$status_log" "retry_timeout" "$trait1" "$trait2" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"
                counter_add "$RETRY_COUNT_FILE" 1
                print_progress_line "RETRY_TIMEOUT" "$trait1" "$trait2" "$runtime" "$attempt"
                attempt=$(( attempt + 1 ))
                sleep 1
                continue
            else
                final_status="timeout"
                write_status_file "$status_log" "$final_status" "$trait1" "$trait2" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"

                counter_add "$TIMEOUT_COUNT_FILE" 1
                counter_add "$DONE_COUNT_FILE" 1
                print_progress_line "TIMEOUT" "$trait1" "$trait2" "$runtime" "$attempt"
                return 0
            fi
        fi

        final_status="failed"
        write_status_file "$status_log" "$final_status" "$trait1" "$trait2" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"

        counter_add "$FAIL_COUNT_FILE" 1
        counter_add "$DONE_COUNT_FILE" 1
        print_progress_line "FAIL" "$trait1" "$trait2" "$runtime" "$attempt"
        return 0
    done
}

print_header() {
    echo "========================================"
    echo "LDSC all-vs-all parallel run"
    echo "Base dir:            $BASE_DIR"
    echo "Munged dir:          $MUNGED_DIR"
    echo "Results dir:         $OUT_DIR"
    echo "Logs dir:            $LOG_DIR"
    echo "Status dir:          $STATUS_DIR"
    echo "Tmp dir:             $TMP_DIR"
    echo "Traits:              $N"
    echo "Pairs:               $TOTAL_PAIRS"
    echo "Parallel jobs:       $JOBS"
    echo "Threads per job:     $THREADS_PER_JOB"
    echo "Timeout per pair:    ${TIMEOUT_SECONDS}s"
    echo "Timeout retries:     $MAX_TIMEOUT_RETRIES"
    echo "Start time:          $(date)"
    echo "========================================"
    echo
}

print_summary() {
    local ok_count fail_count skip_count timeout_count retry_count done_count total_time

    ok_count=$(get_counter "$OK_COUNT_FILE")
    fail_count=$(get_counter "$FAIL_COUNT_FILE")
    skip_count=$(get_counter "$SKIP_COUNT_FILE")
    timeout_count=$(get_counter "$TIMEOUT_COUNT_FILE")
    retry_count=$(get_counter "$RETRY_COUNT_FILE")
    done_count=$(get_counter "$DONE_COUNT_FILE")
    total_time=$(( $(date +%s) - START_ALL ))

    echo
    echo "========================================"
    echo "RUN FINISHED"
    echo "Traits:              $N"
    echo "Pairs:               $TOTAL_PAIRS"
    echo "Done:                $done_count"
    echo "OK:                  $ok_count"
    echo "Failed:              $fail_count"
    echo "Timeout:             $timeout_count"
    echo "Skipped:             $skip_count"
    echo "Timeout retries:     $retry_count"
    echo "Total time:          $(format_seconds "$total_time")"
    echo "End time:            $(date)"
    echo "Results dir:         $OUT_DIR"
    echo "Logs dir:            $LOG_DIR"
    echo "Status dir:          $STATUS_DIR"
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
