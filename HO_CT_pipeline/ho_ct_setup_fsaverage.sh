#!/usr/bin/env bash
# ho_ct_setup_fsaverage.sh
# One-time setup: projects Harvard-Oxford cortical atlas (thr0 and thr25) from
# MNI152 volumetric space onto the fsaverage surface and builds .annot files.
# Outputs are stored in GLOBAL_DIR and reused by all subjects.
# No subject-specific paths are needed.
#
# Usage: bash ho_ct_setup_fsaverage.sh
#
# Load modules before running:
#   module load freesurfer/6.0 fsl

set -euo pipefail

# --- Globals -----------------------------------------------------------------
GLOBAL_DIR="/projectnb/skiran/Computational_R01/2.Scripts/HO_CT_pipeline"
CTAB="${GLOBAL_DIR}/HO48_1based.lut"

[ -f "$CTAB" ] || { echo "ERROR: Missing color table: $CTAB"; exit 1; }

# --- Environment checks ------------------------------------------------------
echo "[INFO] FREESURFER_HOME=${FREESURFER_HOME:-<unset>}"
echo "[INFO] FSLDIR=${FSLDIR:-<unset>}"
[ -n "${FREESURFER_HOME:-}" ] || { echo "ERROR: FREESURFER_HOME not set. Load freesurfer module."; exit 1; }
[ -n "${FSLDIR:-}"          ] || { echo "ERROR: FSLDIR not set. Load fsl module.";                exit 1; }

command -v mri_vol2surf     >/dev/null || { echo "ERROR: mri_vol2surf not found.";     exit 1; }
command -v mri_vol2label    >/dev/null || { echo "ERROR: mri_vol2label not found.";    exit 1; }
command -v mris_label2annot >/dev/null || { echo "ERROR: mris_label2annot not found."; exit 1; }

# mri_vol2surf --mni152reg requires fsaverage in SUBJECTS_DIR.
# mris_label2annot writes output to $SUBJECTS_DIR/fsaverage/label/ with no
# redirect option, so the FreeSurfer installation dir (read-only on SCC) cannot
# be used directly. Copy fsaverage to a writable location under GLOBAL_DIR.
FS_SUBJECTS_READONLY="${FREESURFER_HOME}/subjects"
[ -d "${FS_SUBJECTS_READONLY}/fsaverage/surf" ] || {
  echo "ERROR: fsaverage not found at ${FS_SUBJECTS_READONLY}/fsaverage."
  echo "       Unexpected for a standard FreeSurfer 6.0 installation."
  exit 1
}

ATLAS_DIR="${GLOBAL_DIR}/atlas"
mkdir -p "$ATLAS_DIR"

WRITABLE_SUBJECTS="${ATLAS_DIR}/subjects_tmp"
if [ ! -d "${WRITABLE_SUBJECTS}/fsaverage/surf" ]; then
  echo "[INFO] Copying fsaverage to writable location: ${WRITABLE_SUBJECTS}/fsaverage"
  mkdir -p "$WRITABLE_SUBJECTS"
  cp -r "${FS_SUBJECTS_READONLY}/fsaverage" "${WRITABLE_SUBJECTS}/fsaverage"
else
  echo "[SKIP] Writable fsaverage already present: ${WRITABLE_SUBJECTS}/fsaverage"
fi

export SUBJECTS_DIR="$WRITABLE_SUBJECTS"
echo "[INFO] SUBJECTS_DIR: $SUBJECTS_DIR"

HO_DIR="${FSLDIR}/data/atlases/HarvardOxford"
HO_THR0="${HO_DIR}/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz"
HO_THR25="${HO_DIR}/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz"
[ -f "$HO_THR0"  ] || { echo "ERROR: Missing atlas: $HO_THR0";  exit 1; }
[ -f "$HO_THR25" ] || { echo "ERROR: Missing atlas: $HO_THR25"; exit 1; }

LABEL_DIR="${ATLAS_DIR}/fsaverage_labels"
mkdir -p "${LABEL_DIR}/thr0/lh"  "${LABEL_DIR}/thr0/rh"
mkdir -p "${LABEL_DIR}/thr25/lh" "${LABEL_DIR}/thr25/rh"

# --- Step 1: Project atlas volumes onto fsaverage surface --------------------
# --mni152reg uses $FREESURFER_HOME/average/mni152.register.dat, an affine
# transform between MNI152 and fsaverage (MNI305) space. Appropriate here
# because we are mapping one template to another, not to a patient brain.
# --interp nearest preserves discrete integer label values.
# --projfrac 0.5 samples at the midpoint of the cortical ribbon.

for THR in thr0 thr25; do
  if   [ "$THR" = "thr0"  ]; then ATLAS="$HO_THR0"
  elif [ "$THR" = "thr25" ]; then ATLAS="$HO_THR25"
  fi

  for HEMI in lh rh; do
    SURF_OUT="${ATLAS_DIR}/ho_cort_${THR}_fsaverage_${HEMI}.mgh"
    if [ -f "$SURF_OUT" ]; then
      echo "[SKIP] Already exists: $SURF_OUT"
    else
      echo "[RUN] mri_vol2surf: ${THR} ${HEMI}"
      mri_vol2surf \
        --mov      "$ATLAS" \
        --mni152reg \
        --hemi     "$HEMI" \
        --surf     white \
        --interp   nearest \
        --projfrac 0.5 \
        --o        "$SURF_OUT"
    fi
  done
done

# --- Step 2: Extract per-label surface label files ---------------------------
# mri_vol2label extracts all vertices assigned to integer value IDX from the
# projected surface map into a FreeSurfer .label file.

for THR in thr0 thr25; do
  for HEMI in lh rh; do
    SURF_MAP="${ATLAS_DIR}/ho_cort_${THR}_fsaverage_${HEMI}.mgh"
    LDIR="${LABEL_DIR}/${THR}/${HEMI}"

    ALL_PRESENT=1
    for IDX in $(seq 1 48); do
      [ -f "${LDIR}/${HEMI}.HO_${IDX}.label" ] || { ALL_PRESENT=0; break; }
    done
    if [ "$ALL_PRESENT" = "1" ]; then
      echo "[SKIP] All labels present for ${THR} ${HEMI}"
      continue
    fi

    echo "[RUN] Extracting 48 labels: ${THR} ${HEMI}"
    for IDX in $(seq 1 48); do
      LFILE="${LDIR}/${HEMI}.HO_${IDX}.label"
      [ -f "$LFILE" ] && continue
      mri_vol2label \
        --c    "$SURF_MAP" \
        --id   "$IDX" \
        --surf fsaverage "$HEMI" \
        --l    "$LFILE"
    done
  done
done

# --- Step 3: Assemble per-hemisphere .annot files ----------------------------
# mris_label2annot combines per-label files into a single .annot file using
# the ctab for color and name assignments.
# It writes to $SUBJECTS_DIR/fsaverage/label/ by default; copy to GLOBAL_DIR.
# Vertices not covered by any label are assigned index 0 ("unknown").

for THR in thr0 thr25; do
  for HEMI in lh rh; do
    ANNOT_OUT="${ATLAS_DIR}/${HEMI}.ho_${THR}.annot"
    if [ -f "$ANNOT_OUT" ]; then
      echo "[SKIP] Already exists: $ANNOT_OUT"
      continue
    fi

    echo "[RUN] mris_label2annot: ${THR} ${HEMI}"
    LDIR="${LABEL_DIR}/${THR}/${HEMI}"

    LABEL_ARGS=""
    for IDX in $(seq 1 48); do
      LABEL_ARGS="$LABEL_ARGS --l ${LDIR}/${HEMI}.HO_${IDX}.label"
    done

    # shellcheck disable=SC2086
    mris_label2annot \
      --sd   "$SUBJECTS_DIR" \
      --s    fsaverage \
      --ctab "$CTAB" \
      $LABEL_ARGS \
      --h    "$HEMI" \
      --a    "ho_${THR}"

    SRC="${SUBJECTS_DIR}/fsaverage/label/${HEMI}.ho_${THR}.annot"
    [ -f "$SRC" ] || { echo "ERROR: Expected output not found: $SRC"; exit 1; }
    cp "$SRC" "$ANNOT_OUT"
    echo "[INFO] Copied to: $ANNOT_OUT"
  done
done

# --- Summary -----------------------------------------------------------------
echo
echo "[DONE] Setup complete. Atlas files in: $ATLAS_DIR"
echo
echo "  Surface maps (per-vertex label index on fsaverage):"
for THR in thr0 thr25; do for HEMI in lh rh; do
  echo "    ${ATLAS_DIR}/ho_cort_${THR}_fsaverage_${HEMI}.mgh"
done; done
echo
echo "  Annotation files:"
for THR in thr0 thr25; do for HEMI in lh rh; do
  echo "    ${ATLAS_DIR}/${HEMI}.ho_${THR}.annot"
done; done
echo
echo "QC - verify thr25 annotation on fsaverage in freeview (run as one line):"
echo "  freeview -f \${FREESURFER_HOME}/subjects/fsaverage/surf/lh.inflated:annot=${ATLAS_DIR}/lh.ho_thr25.annot:annot_outline=1"
