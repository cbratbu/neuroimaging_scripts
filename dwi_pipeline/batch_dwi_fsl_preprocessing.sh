#!/usr/bin/env bash
set -euo pipefail

# batch_dwi_fsl_preprocessing.sh
#
# Runs DWI preprocessing on all dwi/ directories found under a cohort
# parent directory. Handles both session and non-session subject structures:
#   <PARENT_DIR>/sub-XXXX/2.preprocessing/dwi/
#   <PARENT_DIR>/sub-XXXX/2.preprocessing/ses-XXXX/dwi/
#
# By default, runs the AP/PA pipeline on subjects with both directions.
# With --synb0, also runs the synb0 pipeline on subjects missing the PA direction.
#
# Usage:
#   batch_dwi_fsl_preprocessing.sh /path/to/cohort_directory [--synb0]

FSL_SCRIPT="/projectnb/openneuro-aphasia/ML_aphasia_scripts/dwi_pipeline/dwi_fsl_preprocessing.sh"
SYNB0_SCRIPT="/projectnb/openneuro-aphasia/ML_aphasia_scripts/dwi_pipeline/dwi_fsl_preprocessing_synb0.sh"

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: $0 /path/to/cohort_directory [--synb0]" >&2
  exit 1
fi

PARENT_DIR="$1"
RUN_SYNB0=false

if [ "${2:-}" = "--synb0" ]; then
  RUN_SYNB0=true
fi

if [ ! -d "$PARENT_DIR" ]; then
  echo "ERROR: Parent directory not found: $PARENT_DIR" >&2
  exit 1
fi

echo "=== DWI FSL preprocessing batch ==="
echo "Parent dir  : ${PARENT_DIR}"
echo "AP/PA script: ${FSL_SCRIPT}"
echo "synb0 script: ${SYNB0_SCRIPT}"
echo "synb0 mode  : ${RUN_SYNB0}"
echo

module purge
module load fsl || true

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
skip_no_pa=0
i=0

for DWI_DIR in "${DWI_DIRS[@]}"; do
  i=$((i+1))

  subj=$(echo "$DWI_DIR" | grep -oP 'sub-[^/]+' | head -n1)
  ses=$(echo "$DWI_DIR" | grep -oP 'ses-[^/]+' | head -n1 || true)
  label="${subj}${ses:+/}${ses:-}"

  eddy_out="${DWI_DIR}/eddy_unwarped_images.nii.gz"

  echo "=== [${i}/${TOTAL}] ${label} ==="
  echo "DWI_DIR : ${DWI_DIR}"

  if [ -f "$eddy_out" ]; then
    echo "[SKIP] eddy_unwarped_images.nii.gz already exists."
    echo
    skip=$((skip+1))
    continue
  fi

  # Detect AP and PA NIfTI presence
  has_ap=false
  has_pa=false
  ls "$DWI_DIR"/*dir-AP_dwi.nii* >/dev/null 2>&1 && has_ap=true
  ls "$DWI_DIR"/*dir-PA_dwi.nii* >/dev/null 2>&1 && has_pa=true

  if $has_ap && $has_pa; then
    echo "AP + PA found. Running AP/PA pipeline..."
    if "$FSL_SCRIPT" "$DWI_DIR"; then
      echo "[OK] ${label}"
      ok=$((ok+1))
    else
      echo "[FAIL] ${label}"
      fail=$((fail+1))
    fi
  elif $has_ap && ! $has_pa; then
    if $RUN_SYNB0; then
      echo "AP only (no PA). Running synb0 pipeline..."
      if "$SYNB0_SCRIPT" "$DWI_DIR"; then
        echo "[OK] ${label}"
        ok=$((ok+1))
      else
        echo "[FAIL] ${label}"
        fail=$((fail+1))
      fi
    else
      echo "[SKIP] AP only — no PA found. Re-run with --synb0 to process these subjects."
      skip_no_pa=$((skip_no_pa+1))
    fi
  else
    echo "[FAIL] No AP NIfTI found in $DWI_DIR"
    fail=$((fail+1))
  fi

  echo
done

echo "=== FSL DWI batch summary ==="
echo "Total subjects  : $TOTAL"
echo "OK              : $ok"
echo "SKIP (done)     : $skip"
echo "SKIP (AP only)  : $skip_no_pa"
echo "FAIL            : $fail"

(( fail == 0 )) || exit 1
