#!/bin/bash
#
# dwi_fsl_preprocessing_synb0.sh
#
# Single-subject DWI preprocessing using FSL + synb0-DisCo:
#   - Builds acqparams.txt and index.txt from AP DWI JSON metadata
#   - Runs synb0-DisCo to synthesize an undistorted b0 from T1w and run TOPUP
#   - Extracts distortion-corrected b0 (hifi_nodif) from TOPUP output and creates
#     brain mask with BET
#   - Runs EDDY (prefers GPU eddy_cuda8.0 if available, otherwise CPU eddy)
#
# For use with single-PE datasets (AP only) that lack a reverse-PE acquisition.
# synb0-DisCo synthesizes an undistorted b0 from the T1w image to use as the
# anatomical target for TOPUP in place of a real PA acquisition. hifi_nodif is
# derived from the distortion-corrected real AP b0 (volume 0 of b0_all_topup).
#
# EXPECTED INPUT DIRECTORY (argument $1):
#   /path/to/sub-XXXXXXX[/ses-YYYY]/dwi
#
# EXPECTED FILES:
#   dwi/  sub-XXXXXXX[_ses-YYYY]_dir-AP_dwi.nii[.gz]
#         sub-XXXXXXX[_ses-YYYY]_dir-AP_dwi.bval
#         sub-XXXXXXX[_ses-YYYY]_dir-AP_dwi.bvec
#         sub-XXXXXXX[_ses-YYYY]_dir-AP_dwi.json  (or pre-existing acqparams.txt)
#   anat/ sub-XXXXXXX[_ses-YYYY]_T1w.nii[.gz]  (sibling of dwi/)
#
# OUTPUTS (in dwi/ unless noted):
#   acqparams.txt, index.txt
#   synb0/INPUTS/  - synb0-DisCo inputs
#   synb0/OUTPUTS/ - synb0-DisCo + TOPUP outputs
#   hifi_nodif.nii.gz, hifi_nodif_brain.nii.gz, hifi_nodif_brain_mask.nii.gz
#   eddy_unwarped_images.nii.gz and associated eddy_* text files
#
# USAGE:
#   ./dwi_fsl_preprocessing_synb0.sh /path/to/sub-XXXXXXX[/ses-YYYY]/dwi
#
# ADJUSTABLE FIELDS:
SYNB0_SIF=/projectnb/openneuro-aphasia/ML_aphasia_scripts/dwi_pipeline/synb0-disco_v3.1.sif
FS_LICENSE=/share/pkg.7/freesurfer/6.0/install/license.txt
#

set -euo pipefail

# -----------------------------
# Argument + directory checks
# -----------------------------
if [ $# -ne 1 ]; then
  echo "Usage: $0 /path/to/sub-XXXXXXX[/ses-YYYY]/dwi"
  exit 1
fi

DWI_DIR="$1"

if [ ! -d "$DWI_DIR" ]; then
  echo "ERROR: Directory not found: $DWI_DIR"
  exit 1
fi

if [ ! -f "$SYNB0_SIF" ]; then
  echo "ERROR: synb0-DisCo SIF not found: $SYNB0_SIF" >&2
  exit 1
fi

if [ ! -f "$FS_LICENSE" ]; then
  echo "ERROR: FreeSurfer license not found: $FS_LICENSE" >&2
  exit 1
fi

ANAT_DIR="$(dirname "$DWI_DIR")/anat"
if [ ! -d "$ANAT_DIR" ]; then
  echo "ERROR: anat/ directory not found at: $ANAT_DIR" >&2
  exit 1
fi

if [ -z "${FSLDIR:-}" ]; then
  echo "ERROR: FSLDIR is not set. Load the FSL module before running this script." >&2
  exit 1
fi

cd "$DWI_DIR"

echo "========================================"
echo " DWI preprocessing (synb0) in: $DWI_DIR"
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
# Detect AP JSON, NIfTI, and T1w
# -----------------------------
AP_JSON=$(ls *dir-AP_dwi.json 2>/dev/null || true)
AP_NII=$(ls *dir-AP_dwi.nii* 2>/dev/null | head -n1 || true)
if [ -z "$AP_NII" ]; then
  echo "ERROR: No dir-AP DWI NIfTI found in $DWI_DIR" >&2
  exit 1
fi
AP_ROOT="${AP_NII%.nii.gz}"
AP_ROOT="${AP_ROOT%.nii}"

T1_NII=$(ls "$ANAT_DIR"/sub-*_T1w.nii 2>/dev/null || ls "$ANAT_DIR"/sub-*_T1w.nii.gz 2>/dev/null || true)
if [ -z "$T1_NII" ]; then
  echo "ERROR: No T1w file found in $ANAT_DIR" >&2
  exit 1
fi

echo "AP NIfTI: $AP_NII"
echo "T1w:      $T1_NII"
echo ""

# -----------------------------
# Build acqparams.txt
# Uses existing acqparams.txt if present, otherwise builds from JSON.
# Row 1: real AP b0 with actual TotalReadoutTime
# Row 2: synthetic undistorted b0 - same PE direction, TRO=0
# -----------------------------
if [ -f acqparams.txt ]; then
  echo "Using existing acqparams.txt:"
  cat acqparams.txt
  echo ""
else
  echo "Building acqparams.txt ..."
  if [ -z "$AP_JSON" ]; then
    echo "ERROR: No dir-AP JSON found and no acqparams.txt present. Cannot build acqparams.txt." >&2
    exit 1
  fi
  AP_PE=$(grep -o '"PhaseEncodingDirection": *"[^"]*"' "$AP_JSON" | awk -F'"' '{print $4}')
  AP_TRO=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$AP_JSON" | awk '{print $2}')
  AP_VEC=$(pevec "$AP_PE")
  echo "$AP_VEC $AP_TRO" >  acqparams.txt
  echo "$AP_VEC 0"        >> acqparams.txt
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
# Step 1: Prepare synb0 inputs
# -----------------------------
echo "Step 1: Preparing synb0-DisCo inputs ..."

mkdir -p synb0/INPUTS synb0/OUTPUTS

fslroi "$AP_NII" synb0/INPUTS/b0.nii.gz 0 1
fslchfiletype NIFTI_GZ "$T1_NII" synb0/INPUTS/T1.nii.gz

cp acqparams.txt synb0/INPUTS/acqparams.txt

echo "synb0 INPUTS:"
ls -1 synb0/INPUTS/
echo ""

# -----------------------------
# Step 2: Run synb0-DisCo (synthesis + TOPUP)
# -----------------------------
echo "Step 2: synb0-DisCo ..."

singularity run -e \
  -B "$(pwd)/synb0/INPUTS":/INPUTS \
  -B "$(pwd)/synb0/OUTPUTS":/OUTPUTS \
  -B "$FS_LICENSE":/extra/freesurfer/license.txt \
  "$SYNB0_SIF"

echo "synb0-DisCo finished. OUTPUTS:"
ls -1 synb0/OUTPUTS/
echo ""

if [ ! -f synb0/OUTPUTS/topup_fieldcoef.nii.gz ]; then
  echo "ERROR: TOPUP did not produce expected outputs in synb0/OUTPUTS/" >&2
  exit 1
fi

# -----------------------------
# Step 3: Extract hifi_nodif + BET
# hifi_nodif is taken from volume 0 of b0_all_topup, which is the
# distortion-corrected real AP b0. Volume 1 is the corrected synthetic b0
# and is not averaged in, as mixing real and synthesized signal is not appropriate.
# -----------------------------
echo "Step 3: hifi_nodif + BET ..."

fslroi synb0/OUTPUTS/b0_all_topup hifi_nodif 0 1
bet hifi_nodif hifi_nodif_brain -m -f 0.2

echo "Generated hifi_nodif/BET outputs:"
ls -1 hifi_nodif*
echo ""

# -----------------------------
# Step 4: Choose EDDY backend
# -----------------------------
echo "Step 4: EDDY ..."

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
  --topup=synb0/OUTPUTS/topup \
  --flm=quadratic \
  --out=eddy_unwarped_images \
  --data_is_shelled

echo ""
echo "EDDY finished. Outputs:"
ls -1 eddy_unwarped_images*
echo "========================================"
echo " DWI preprocessing completed for: $DWI_DIR"
echo "========================================"