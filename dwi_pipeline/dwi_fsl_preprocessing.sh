#!/bin/bash
#
# dwi_fsl_preprocessing.sh
#
# Single-subject DWI preprocessing using FSL:
#   - Builds acqparams.txt and index.txt from AP/PA JSON metadata
#   - Extracts b0 volumes from AP/PA DWI and merges them
#   - Runs TOPUP to estimate and correct susceptibility distortion
#   - Averages corrected b0 volumes (hifi_nodif)
#   - Runs BET to create a brain mask from hifi_nodif
#   - Runs EDDY (prefers GPU eddy_cuda8.0 if a GPU is available, otherwise CPU eddy)
#
# EXPECTED INPUT DIRECTORY (argument $1):
#   /path/to/sub-XXXXXXX/dwi
#
# EXPECTED FILES IN THAT DIRECTORY:
#   sub-XXXXXXX_dir-AP_dwi.nii(.gz)
#   sub-XXXXXXX_dir-AP_dwi.bval
#   sub-XXXXXXX_dir-AP_dwi.bvec
#   sub-XXXXXXX_dir-AP_dwi.json  (or pre-existing acqparams.txt)
#   sub-XXXXXXX_dir-PA_dwi.nii(.gz)
#   sub-XXXXXXX_dir-PA_dwi.bval
#   sub-XXXXXXX_dir-PA_dwi.bvec
#   sub-XXXXXXX_dir-PA_dwi.json  (or pre-existing acqparams.txt)
#
# OUTPUTS (in the same directory):
#   acqparams.txt, index.txt
#   AP.nii.gz, PA.nii.gz, AP_PA_b0.nii.gz
#   topup_AP_PA_b0_* (fieldcoef, fout, iout, movpar)
#   hifi_nodif.nii.gz, hifi_nodif_brain.nii.gz, hifi_nodif_brain_mask.nii.gz
#   eddy_unwarped_images.nii.gz and associated eddy_* text files
#
# USAGE:
#   ./dwi_fsl_preprocessing.sh /path/to/sub-XXXXXXX/dwi
#

set -euo pipefail

# -----------------------------
# Argument + directory checks
# -----------------------------
if [ $# -ne 1 ]; then
  echo "Usage: $0 /path/to/sub-XXXXXXX/dwi"
  exit 1
fi

DWI_DIR="$1"

if [ ! -d "$DWI_DIR" ]; then
  echo "ERROR: Directory not found: $DWI_DIR"
  exit 1
fi

if [ -z "${FSLDIR:-}" ]; then
  echo "ERROR: FSLDIR is not set. Load the FSL module before running this script." >&2
  exit 1
fi

cd "$DWI_DIR"

echo "========================================"
echo " DWI preprocessing in: $DWI_DIR"
echo "========================================"

# -----------------------------
# Helper: map PhaseEncodingDirection -> vector
# -----------------------------
pevec() {
  local phase="$1"
  case "$phase" in
    i)  echo "1 0 0" ;;
    i-) echo "-1 0 0" ;;
    j)  echo "0 1 0" ;;
    j-) echo "0 -1 0" ;;
    k)  echo "0 0 1" ;;
    k-) echo "0 0 -1" ;;
    *)  echo "ERROR: Unknown PhaseEncodingDirection: $phase" >&2; exit 1 ;;
  esac
}

# -----------------------------
# Detect AP/PA JSON & NIfTI roots
# -----------------------------
AP_JSON=$(ls *dir-AP_dwi.json 2>/dev/null || true)
PA_JSON=$(ls *dir-PA_dwi.json 2>/dev/null || true)

AP_NII=$(ls *dir-AP_dwi.nii* 2>/dev/null | head -n1 || true)
PA_NII=$(ls *dir-PA_dwi.nii* 2>/dev/null | head -n1 || true)

if [ -z "$AP_NII" ]; then
  echo "ERROR: No dir-AP DWI NIfTI found in $DWI_DIR" >&2
  exit 1
fi
if [ -z "$PA_NII" ]; then
  echo "ERROR: No dir-PA DWI NIfTI found in $DWI_DIR" >&2
  exit 1
fi

AP_ROOT="${AP_NII%.nii.gz}"
AP_ROOT="${AP_ROOT%.nii}"
PA_ROOT="${PA_NII%.nii.gz}"
PA_ROOT="${PA_ROOT%.nii}"

echo "AP NIfTI: $AP_NII"
echo "PA NIfTI: $PA_NII"
echo ""

# -----------------------------
# Build acqparams.txt
# Uses existing acqparams.txt if present, otherwise builds from JSONs.
# -----------------------------
if [ -f acqparams.txt ]; then
  echo "Using existing acqparams.txt:"
  cat acqparams.txt
  echo ""
else
  echo "Building acqparams.txt ..."
  if [ -z "$AP_JSON" ] || [ -z "$PA_JSON" ]; then
    echo "ERROR: No acqparams.txt found and AP or PA JSON missing. Cannot build acqparams.txt." >&2
    exit 1
  fi
  AP_PE=$(grep -o '"PhaseEncodingDirection": *"[^"]*"' "$AP_JSON" | awk -F'"' '{print $4}')
  AP_TRO=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$AP_JSON" | awk '{print $2}')
  PA_PE=$(grep -o '"PhaseEncodingDirection": *"[^"]*"' "$PA_JSON" | awk -F'"' '{print $4}')
  PA_TRO=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$PA_JSON" | awk '{print $2}')
  AP_VEC=$(pevec "$AP_PE")
  PA_VEC=$(pevec "$PA_PE")
  echo "$AP_VEC $AP_TRO" >  acqparams.txt
  echo "$PA_VEC $PA_TRO" >> acqparams.txt
  echo "acqparams.txt:"
  cat acqparams.txt
  echo ""
fi

# -----------------------------
# Build index.txt
# -----------------------------
echo "Building index.txt ..."

NVOLS=$(fslval "$AP_NII" dim4)
if [ -z "$NVOLS" ]; then
  echo "ERROR: Could not read number of volumes from $AP_NII" >&2
  exit 1
fi

{
  for ((i=1; i<=NVOLS; i++)); do
    echo -n "1 "
  done
  echo ""
} > index.txt

echo "index.txt (first 80 characters):"
head -c 80 index.txt; echo ""
echo "Total volumes (AP dim4): $NVOLS"
echo ""

# -----------------------------
# Step 1: Extract b0 and merge
# -----------------------------
echo "Step 1: fslroi + fslmerge ..."

fslroi "$AP_NII" AP 0 1
fslroi "$PA_NII" PA 0 1
fslmerge -t AP_PA_b0 AP PA

echo "AP dim4:       $(fslval AP dim4)"
echo "PA dim4:       $(fslval PA dim4)"
echo "AP_PA_b0 dim4: $(fslval AP_PA_b0 dim4)"
echo ""

# -----------------------------
# Step 2: TOPUP + mean b0 + BET
# -----------------------------
echo "Step 2: TOPUP + hifi_nodif + BET ..."

TOPUP_CFG="$FSLDIR/etc/flirtsch/b02b0.cnf"
if [ ! -f "$TOPUP_CFG" ]; then
  echo "ERROR: Cannot find b02b0.cnf at $TOPUP_CFG" >&2
  exit 1
fi

topup \
  --imain=AP_PA_b0 \
  --datain=acqparams.txt \
  --config="$TOPUP_CFG" \
  --out=topup_AP_PA_b0 \
  --iout=topup_AP_PA_b0_iout \
  --fout=topup_AP_PA_b0_fout

fslmaths topup_AP_PA_b0_iout -Tmean hifi_nodif
bet hifi_nodif hifi_nodif_brain -m -f 0.2

echo "Generated TOPUP/BET outputs:"
ls -1 topup_AP_PA_b0* hifi_nodif*
echo ""

# -----------------------------
# Step 3: Choose EDDY backend
# -----------------------------
echo "Step 3: EDDY ..."

EDDY_CMD="eddy"

if command -v nvidia-smi >/dev/null 2>&1; then
  if nvidia-smi -L >/dev/null 2>&1; then
    module load cuda/8.0 >/dev/null 2>&1 || true
    if [ -x "$FSLDIR/bin/eddy_cuda8.0" ]; then
      EDDY_CMD="eddy_cuda8.0"
    fi
  fi
fi

echo "Using EDDY command: $EDDY_CMD"
echo ""

# -----------------------------
# Run EDDY
# -----------------------------
$EDDY_CMD \
  --imain="$AP_ROOT" \
  --mask=hifi_nodif_brain_mask \
  --index=index.txt \
  --acqp=acqparams.txt \
  --bvecs="${AP_ROOT}.bvec" \
  --bvals="${AP_ROOT}.bval" \
  --fwhm=0 \
  --topup=topup_AP_PA_b0 \
  --flm=quadratic \
  --out=eddy_unwarped_images \
  --data_is_shelled

echo ""
echo "EDDY finished. Outputs:"
ls -1 eddy_unwarped_images*
echo "========================================"
echo " DWI preprocessing completed for: $DWI_DIR"
echo "========================================"