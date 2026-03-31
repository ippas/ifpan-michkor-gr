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

H2_BASE_DIR="${BASE_DIR}/results_h2"

OUT_DIR="${H2_BASE_DIR}/results"
LOG_DIR="${H2_BASE_DIR}/logs"
STATUS_DIR="${H2_BASE_DIR}/status"
TMP_DIR="${H2_BASE_DIR}/tmp"

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
    echo "  bash run_ldsc_h2_all_parallel_v1.sh [options]"
    echo
    echo "Options:"
    echo "  -j, --jobs N                Number of parallel jobs (default: 10)"
    echo "  -t, --timeout-seconds N     Timeout per trait in seconds (default: 300)"
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

START_ALL=0
N=0

LOCK_FILE=""
DONE_COUNT_FILE=""
OK_COUNT_FILE=""
FAIL_COUNT_FILE=""
SKIP_COUNT_FILE=""
TIMEOUT_COUNT_FILE=""
RETRY_COUNT_FILE=""

FILES=()

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
    mkdir -p "$H2_BASE_DIR"
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

    if [[ "$N" -lt 1 ]]; then
        echo "ERROR: no munged files found"
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
    local trait="$2"
    local runtime="$3"
    local attempt="$4"

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

        left_count=$(( N - done_count ))
        percent=$(awk -v d="$done_count" -v t="$N" 'BEGIN {printf "%.2f", 100*d/t}')

        if [[ "$runtime" -ge 0 ]]; then
            echo "[${done_count}/${N} | ${percent}% | left:${left_count} | ok:${ok_count} fail:${fail_count} timeout:${timeout_count} skip:${skip_count} retry:${retry_count}] ${label} ${trait} (attempt:${attempt}, $(format_seconds "$runtime"))"
        else
            echo "[${done_count}/${N} | ${percent}% | left:${left_count} | ok:${ok_count} fail:${fail_count} timeout:${timeout_count} skip:${skip_count} retry:${retry_count}] ${label} ${trait} (attempt:${attempt})"
        fi
    ) 200>"$LOCK_FILE"
}

write_status_file() {
    local status_log="$1"
    local status="$2"
    local trait="$3"
    local runtime="$4"
    local exit_code="$5"
    local out_prefix="$6"
    local stdout_log="$7"
    local stderr_log="$8"
    local attempt="$9"

    {
        echo "status=${status}"
        echo "trait=${trait}"
        echo "runtime_sec=${runtime}"
        echo "exit_code=${exit_code}"
        echo "output_prefix=${out_prefix}"
        echo "stdout_log=${stdout_log}"
        echo "stderr_log=${stderr_log}"
        echo "attempt=${attempt}"
        echo "timestamp=$(date '+%F %T')"
    } > "$status_log"
}

run_trait() {
    local file="$1"

    local trait
    local out_prefix status_log done_marker
    local start_one end_one runtime
    local attempt max_attempts
    local stdout_log stderr_log
    local exit_code final_status

    trait=$(basename "$file" .sumstats.gz)

    out_prefix="${OUT_DIR}/${trait}"
    status_log="${STATUS_DIR}/${trait}.status.txt"
    done_marker="${STATUS_DIR}/${trait}.done"

    if [[ -f "$done_marker" ]] || { [[ -f "$status_log" ]] && grep -q '^status=ok$' "$status_log"; }; then
        counter_add "$SKIP_COUNT_FILE" 1
        counter_add "$DONE_COUNT_FILE" 1
        print_progress_line "SKIP" "$trait" 0 0
        return 0
    fi

    max_attempts=$(( MAX_TIMEOUT_RETRIES + 1 ))
    attempt=1
    final_status="failed"

    while [[ "$attempt" -le "$max_attempts" ]]; do
        stdout_log="${LOG_DIR}/${trait}.attempt${attempt}.out.log"
        stderr_log="${LOG_DIR}/${trait}.attempt${attempt}.err.log"

        start_one=$(date +%s)

        OMP_NUM_THREADS="$THREADS_PER_JOB" \
        OPENBLAS_NUM_THREADS="$THREADS_PER_JOB" \
        MKL_NUM_THREADS="$THREADS_PER_JOB" \
        NUMEXPR_NUM_THREADS="$THREADS_PER_JOB" \
        VECLIB_MAXIMUM_THREADS="$THREADS_PER_JOB" \
        timeout --signal=TERM --kill-after=30s "${TIMEOUT_SECONDS}s" \
        python "$LDSC" \
            --h2 "$file" \
            --ref-ld-chr "$REF" \
            --w-ld-chr "$WLD" \
            --out "$out_prefix" \
            > "$stdout_log" 2> "$stderr_log"

        exit_code=$?

        end_one=$(date +%s)
        runtime=$(( end_one - start_one ))

        if [[ "$exit_code" -eq 0 ]]; then
            final_status="ok"
            write_status_file "$status_log" "$final_status" "$trait" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"
            touch "$done_marker"

            counter_add "$OK_COUNT_FILE" 1
            counter_add "$DONE_COUNT_FILE" 1
            print_progress_line "OK  " "$trait" "$runtime" "$attempt"
            return 0
        fi

        if [[ "$exit_code" -eq 124 || "$exit_code" -eq 137 ]]; then
            if [[ "$attempt" -lt "$max_attempts" ]]; then
                write_status_file "$status_log" "retry_timeout" "$trait" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"
                counter_add "$RETRY_COUNT_FILE" 1
                print_progress_line "RETRY_TIMEOUT" "$trait" "$runtime" "$attempt"
                attempt=$(( attempt + 1 ))
                sleep 1
                continue
            else
                final_status="timeout"
                write_status_file "$status_log" "$final_status" "$trait" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"

                counter_add "$TIMEOUT_COUNT_FILE" 1
                counter_add "$DONE_COUNT_FILE" 1
                print_progress_line "TIMEOUT" "$trait" "$runtime" "$attempt"
                return 0
            fi
        fi

        final_status="failed"
        write_status_file "$status_log" "$final_status" "$trait" "$runtime" "$exit_code" "$out_prefix" "$stdout_log" "$stderr_log" "$attempt"

        counter_add "$FAIL_COUNT_FILE" 1
        counter_add "$DONE_COUNT_FILE" 1
        print_progress_line "FAIL" "$trait" "$runtime" "$attempt"
        return 0
    done
}

print_header() {
    echo "========================================"
    echo "LDSC heritability parallel run"
    echo "Base dir:            $BASE_DIR"
    echo "Munged dir:          $MUNGED_DIR"
    echo "H2 base dir:         $H2_BASE_DIR"
    echo "Results dir:         $OUT_DIR"
    echo "Logs dir:            $LOG_DIR"
    echo "Status dir:          $STATUS_DIR"
    echo "Tmp dir:             $TMP_DIR"
    echo "Traits:              $N"
    echo "Parallel jobs:       $JOBS"
    echo "Threads per job:     $THREADS_PER_JOB"
    echo "Timeout per trait:   ${TIMEOUT_SECONDS}s"
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
    echo "Done:                $done_count"
    echo "OK:                  $ok_count"
    echo "Failed:              $fail_count"
    echo "Timeout:             $timeout_count"
    echo "Skipped:             $skip_count"
    echo "Timeout retries:     $retry_count"
    echo "Total time:          $(format_seconds "$total_time")"
    echo "End time:            $(date)"
    echo "H2 base dir:         $H2_BASE_DIR"
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
print_header

for file in "${FILES[@]}"; do
    wait_for_slot
    run_trait "$file" &
done

wait
print_summary
