#!/usr/bin/env bash
set -euo pipefail

# batch_dwi_dsi_processing.sh
#
# Runs dwi_dsi_pipeline.sh on all dwi/ directories found under a cohort
# parent directory. Handles both session and non-session subject structures:
#   <PARENT_DIR>/sub-XXXX/2.preprocessing/dwi/
#   <PARENT_DIR>/sub-XXXX/2.preprocessing/ses-XXXX/dwi/
#
# Usage:
#   batch_dwi_dsi_processing.sh /path/to/cohort_directory

DSI_SCRIPT="/projectnb/openneuro-aphasia/ML_aphasia_scripts/dwi_pipeline/dwi_dsi_pipeline.sh"

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /path/to/cohort_directory" >&2
  exit 1
fi

PARENT_DIR="$1"

if [ ! -d "$PARENT_DIR" ]; then
  echo "ERROR: Parent directory not found: $PARENT_DIR" >&2
  exit 1
fi

echo "=== DWI DSI Studio batch ==="
echo "Parent dir : ${PARENT_DIR}"
echo "Script     : ${DSI_SCRIPT}"
echo

module purge || true

mapfile -t DWI_DIRS < <(
  {
    find "$PARENT_DIR" -maxdepth 3 -type d -name "dwi" -path "*/2.preprocessing/dwi"
    find "$PARENT_DIR" -maxdepth 4 -type d -name "dwi" -path "*/2.preprocessing/ses-*/dwi"
  } | sort -u
)
TOTAL=${#DWI_DIRS[@]}

echo "Found $TOTAL dwi directories:"
printf '  %s\n' "${DWI_DIRS[@]}"
echo

ok=0
fail=0
skip=0
i=0

for DWI_DIR in "${DWI_DIRS[@]}"; do
  i=$((i+1))

  subj=$(echo "$DWI_DIR" | grep -oP 'sub-[^/]+' | head -n1)
  ses=$(echo "$DWI_DIR" | grep -oP 'ses-[^/]+' | head -n1 || true)
  label="${subj}${ses:+/}${ses:-}"

  eddy_out="${DWI_DIR}/eddy_unwarped_images.nii.gz"

  echo "=== [${i}/${TOTAL}] ${label} ==="
  echo "DWI_DIR : ${DWI_DIR}"

  if [ ! -f "$eddy_out" ]; then
    echo "[SKIP] eddy_unwarped_images.nii.gz not found — FSL preprocessing not yet run."
    echo
    skip=$((skip+1))
    continue
  fi

  bval_file=$(ls "${DWI_DIR}"/sub-*_dir-AP_dwi.bval 2>/dev/null | head -n1 || true)
  if [ -n "$bval_file" ]; then
    subject_id=$(basename "$bval_file" | sed 's/_dir-AP_dwi.bval//')
    final_tsv="${DWI_DIR}/tract_stats/${subject_id}_eddy_unwarped_images_autotrack_stats.tsv"
    if [ -f "$final_tsv" ]; then
      echo "[SKIP] Final TSV already exists: ${final_tsv}"
      echo
      skip=$((skip+1))
      continue
    fi
  fi

  # Build skip flags based on existing outputs
  SKIP_FLAGS=()

  sz_file=$(ls "${DWI_DIR}/${subject_id}_eddy_unwarped_images.sz" 2>/dev/null || true)
  if [ -f "$sz_file" ]; then
    SKIP_FLAGS+=("--skip-src")
    echo "[INFO] SRC exists, will skip SRC generation."
  fi

  fib_file="${DWI_DIR}/${subject_id}_eddy_unwarped_images.fib.gz"
  if [ -f "$fib_file" ]; then
    SKIP_FLAGS+=("--skip-rec")
    echo "[INFO] FIB exists, will skip GQI reconstruction."
  fi

  echo "Running DSI Studio pipeline..."
  if "$DSI_SCRIPT" "${SKIP_FLAGS[@]}" "$DWI_DIR"; then
    echo "[OK] ${label}"
    ok=$((ok+1))
  else
    echo "[FAIL] ${label}"
    fail=$((fail+1))
  fi

  echo
done

echo "=== DSI Studio batch summary ==="
echo "Total : $TOTAL"
echo "OK    : $ok"
echo "SKIP  : $skip"
echo "FAIL  : $fail"

(( fail == 0 )) || exit 1
