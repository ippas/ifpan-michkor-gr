#!/usr/bin/env bash

set -u
set -o pipefail

############################################
# CONFIG
############################################

BASE_DIR="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_analysis/main_sampleSize10000noMissingSubcategory"

MDD="${BASE_DIR}/munged/ieu-a-1188.sumstats.gz"
MUNGED_DIR="${BASE_DIR}/munged"

LDSC="/home/mateusz/projects/ifpan-michkor-gr/tools/ldsc/ldsc.py"
REF="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"
WLD="/home/mateusz/projects/ifpan-michkor-gr/data/ldsc_ref/eur_w_ld_chr/"

############################################
# FUNCTIONS
############################################

format_seconds() {
    local total_seconds=$1

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

############################################
# CHECK INPUTS
############################################

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: BASE_DIR not found:"
    echo "$BASE_DIR"
    exit 1
fi

if [[ ! -f "$MDD" ]]; then
    echo "ERROR: reference file not found:"
    echo "$MDD"
    exit 1
fi

if [[ ! -d "$MUNGED_DIR" ]]; then
    echo "ERROR: munged directory not found:"
    echo "$MUNGED_DIR"
    exit 1
fi

if [[ ! -f "$LDSC" ]]; then
    echo "ERROR: ldsc.py not found:"
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

############################################
# PREP OUTPUT DIRS
############################################

REF_NAME=$(basename "$MDD" .sumstats.gz)

OUT_DIR="${BASE_DIR}/rg_results_${REF_NAME}"
LOG_DIR="${BASE_DIR}/rg_logs_${REF_NAME}"
STATUS_DIR="${BASE_DIR}/rg_status_${REF_NAME}"

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$STATUS_DIR"

############################################
# PREP FILE LIST
############################################

REF_BASENAME=$(basename "$MDD")

FILES=()
while IFS= read -r -d '' file; do
    if [[ "$(basename "$file")" != "$REF_BASENAME" ]]; then
        FILES+=("$file")
    fi
done < <(find "$MUNGED_DIR" -maxdepth 1 -type f -name "*.sumstats.gz" -print0 | sort -z)

TOTAL=${#FILES[@]}

if [[ "$TOTAL" -eq 0 ]]; then
    echo "No comparison files found."
    exit 0
fi

DONE=0
FAILED=0
SKIPPED=0
SUCCESS_TIMED=0
SUM_TIME=0

START_ALL=$(date +%s)

############################################
# HEADER
############################################

echo "========================================"
echo "LDSC one-vs-all"
echo "Reference trait: $REF_NAME"
echo "Reference file:  $MDD"
echo "Munged dir:      $MUNGED_DIR"
echo "Results dir:     $OUT_DIR"
echo "Logs dir:        $LOG_DIR"
echo "Status dir:      $STATUS_DIR"
echo "Total traits:    $TOTAL"
echo "Start time:      $(date)"
echo "========================================"
echo

############################################
# LOOP
############################################

for FILE in "${FILES[@]}"; do
    TRAIT=$(basename "$FILE" .sumstats.gz)

    OUT_PREFIX="${OUT_DIR}/${REF_NAME}_vs_${TRAIT}"
    STDOUT_LOG="${LOG_DIR}/${REF_NAME}_vs_${TRAIT}.out.log"
    STDERR_LOG="${LOG_DIR}/${REF_NAME}_vs_${TRAIT}.err.log"
    STATUS_LOG="${STATUS_DIR}/${REF_NAME}_vs_${TRAIT}.status.txt"

    COMPLETED_TOTAL=$((DONE + SKIPPED))

    if [[ -f "${OUT_PREFIX}.log" ]]; then
        SKIPPED=$((SKIPPED + 1))
        COMPLETED_TOTAL=$((DONE + SKIPPED))
        REMAINING=$((TOTAL - COMPLETED_TOTAL))
        ELAPSED=$(( $(date +%s) - START_ALL ))

        if [[ "$SUCCESS_TIMED" -gt 0 ]]; then
            AVG_TIME=$((SUM_TIME / SUCCESS_TIMED))
            ETA=$((AVG_TIME * REMAINING))
            AVG_FMT=$(format_seconds "$AVG_TIME")
            ETA_FMT=$(format_seconds "$ETA")
        else
            AVG_FMT="NA"
            ETA_FMT="NA"
        fi

        PERCENT=$(awk "BEGIN {printf \"%.2f\", (${COMPLETED_TOTAL}/${TOTAL})*100}")

        echo "[$COMPLETED_TOTAL/$TOTAL] SKIP   $TRAIT"
        echo "  progress: ${PERCENT}%"
        echo "  elapsed:  $(format_seconds "$ELAPSED")"
        echo "  avg time: ${AVG_FMT}"
        echo "  ETA:      ${ETA_FMT}"
        echo "  failed:   $FAILED"
        echo "  skipped:  $SKIPPED"
        echo "----------------------------------------"

        {
            echo "status=skipped"
            echo "trait=${TRAIT}"
            echo "reference=${REF_NAME}"
        } > "$STATUS_LOG"

        continue
    fi

    START_ONE=$(date +%s)

    python "$LDSC" \
        --rg "${MDD},${FILE}" \
        --ref-ld-chr "$REF" \
        --w-ld-chr "$WLD" \
        --out "$OUT_PREFIX" \
        > "$STDOUT_LOG" 2> "$STDERR_LOG"

    EXIT_CODE=$?

    END_ONE=$(date +%s)
    RUNTIME=$((END_ONE - START_ONE))

    DONE=$((DONE + 1))
    COMPLETED_TOTAL=$((DONE + SKIPPED))

    if [[ "$EXIT_CODE" -eq 0 ]]; then
        STATUS="OK"
        SUM_TIME=$((SUM_TIME + RUNTIME))
        SUCCESS_TIMED=$((SUCCESS_TIMED + 1))
    else
        STATUS="FAIL"
        FAILED=$((FAILED + 1))
    fi

    REMAINING=$((TOTAL - COMPLETED_TOTAL))
    ELAPSED=$((END_ONE - START_ALL))

    if [[ "$SUCCESS_TIMED" -gt 0 ]]; then
        AVG_TIME=$((SUM_TIME / SUCCESS_TIMED))
        ETA=$((AVG_TIME * REMAINING))
        AVG_FMT=$(format_seconds "$AVG_TIME")
        ETA_FMT=$(format_seconds "$ETA")
    else
        AVG_FMT="NA"
        ETA_FMT="NA"
    fi

    PERCENT=$(awk "BEGIN {printf \"%.2f\", (${COMPLETED_TOTAL}/${TOTAL})*100}")

    echo "[$COMPLETED_TOTAL/$TOTAL] ${STATUS}    $TRAIT"
    echo "  runtime:  $(format_seconds "$RUNTIME")"
    echo "  progress: ${PERCENT}%"
    echo "  elapsed:  $(format_seconds "$ELAPSED")"
    echo "  avg time: ${AVG_FMT}"
    echo "  ETA:      ${ETA_FMT}"
    echo "  failed:   $FAILED"
    echo "  skipped:  $SKIPPED"
    echo "----------------------------------------"

    {
        echo "status=$(echo "$STATUS" | tr '[:upper:]' '[:lower:]')"
        echo "trait=${TRAIT}"
        echo "reference=${REF_NAME}"
        echo "runtime_sec=${RUNTIME}"
        echo "exit_code=${EXIT_CODE}"
        echo "stdout_log=${STDOUT_LOG}"
        echo "stderr_log=${STDERR_LOG}"
        echo "output_prefix=${OUT_PREFIX}"
    } > "$STATUS_LOG"
done

############################################
# SUMMARY
############################################

END_ALL=$(date +%s)
TOTAL_TIME=$((END_ALL - START_ALL))
COMPLETED_TOTAL=$((DONE + SKIPPED))

echo
echo "========================================"
echo "RUN FINISHED"
echo "Reference trait: $REF_NAME"
echo "Completed total: $COMPLETED_TOTAL / $TOTAL"
echo "Executed:        $DONE"
echo "Failed:          $FAILED"
echo "Skipped:         $SKIPPED"
echo "Total time:      $(format_seconds "$TOTAL_TIME")"
if [[ "$SUCCESS_TIMED" -gt 0 ]]; then
    echo "Mean runtime:    $(format_seconds $((SUM_TIME / SUCCESS_TIMED)))"
else
    echo "Mean runtime:    NA"
fi
echo "Results dir:     $OUT_DIR"
echo "Logs dir:        $LOG_DIR"
echo "Status dir:      $STATUS_DIR"
echo "End time:        $(date)"
echo "========================================"
