#!/usr/bin/env bash
# batch_build_cfgs.sh
set -euo pipefail

EXTRA_FLAGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --useFLAIR|--useT2w)
      EXTRA_FLAGS+=("$1")
      shift
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 [--useFLAIR] [--useT2w] {PWA|HC} /path/to/parent_dir" >&2
  exit 1
fi

POP="$1"
PARENT_DIR="${2%/}"

if [[ "$POP" != "PWA" && "$POP" != "HC" ]]; then
  echo "ERROR: first arg must be PWA or HC (got: $POP)" >&2
  exit 1
fi

if [[ ! -d "$PARENT_DIR" ]]; then
  echo "ERROR: parent dir not found: $PARENT_DIR" >&2
  exit 1
fi

SCRIPT_DIR="/projectnb/openneuro-aphasia/ML_aphasia_scripts/CONN_Preprocessing"

if [[ "$POP" == "PWA" ]]; then
  BUILDER="${SCRIPT_DIR}/build_conn_cfg_PWA_single.sh"
else
  BUILDER="${SCRIPT_DIR}/build_conn_cfg_HC_single.sh"
fi

if [[ ! -x "$BUILDER" ]]; then
  echo "ERROR: builder not found or not executable: $BUILDER" >&2
  exit 1
fi

mapfile -t SUBJECTS < <(find "$PARENT_DIR" -maxdepth 1 -type d -name "sub-*" | sort)
TOTAL_SUBJECTS=${#SUBJECTS[@]}

echo "=== Batch build CONN cfgs (${POP}) ==="
echo "Parent dir : ${PARENT_DIR}"
echo "Builder    : ${BUILDER}"
echo "Extra flags: ${EXTRA_FLAGS[*]}"
echo "Found      : ${TOTAL_SUBJECTS} subjects"
echo

processed_units=0
skipped_units=0
total_units=0

failed_subjects=0
failed_list=()

i=0
for SUBJ_DIR in "${SUBJECTS[@]}"; do
  i=$((i+1))
  subj="$(basename "$SUBJ_DIR")"

  echo "=== [${i}/${TOTAL_SUBJECTS}] ${subj} ==="

  out="$("$BUILDER" "${EXTRA_FLAGS[@]}" "$SUBJ_DIR" 2>&1)" || {
    echo "$out"
    echo "[FAIL] ${subj}: builder exited nonzero"
    failed_subjects=$((failed_subjects+1))
    failed_list+=("$subj")
    echo
    continue
  }

  echo "$out"

  wrote=0
  skipped=0
  line="$(printf "%s\n" "$out" | grep -E '^Done\. Wrote: [0-9]+[[:space:]]+Skipped: [0-9]+' | tail -n 1 || true)"
  if [[ -n "$line" ]]; then
    wrote="$(printf "%s" "$line" | sed -E 's/^Done\. Wrote: ([0-9]+)[[:space:]]+Skipped: ([0-9]+).*$/\1/')"
    skipped="$(printf "%s" "$line" | sed -E 's/^Done\. Wrote: ([0-9]+)[[:space:]]+Skipped: ([0-9]+).*$/\2/')"
  fi

  processed_units=$((processed_units + wrote))
  skipped_units=$((skipped_units + skipped))
  total_units=$((total_units + wrote + skipped))

  echo
done

echo "=== Batch cfg build summary ==="
echo "Processed units : ${processed_units} / ${total_units}"
echo "Skipped units   : ${skipped_units}"
echo "Failed subjects : ${failed_subjects}"
if (( failed_subjects > 0 )); then
  echo "Failed list     : ${failed_list[*]}"
fi

(( failed_subjects == 0 )) || exit 1
