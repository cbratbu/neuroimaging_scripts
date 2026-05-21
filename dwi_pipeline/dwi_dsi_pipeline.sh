#!/usr/bin/env bash
set -euo pipefail

# dwi_dsi_pipeline.sh
#
# DSI Studio pipeline: SRC generation, GQI reconstruction, AutoTrack tractography,
# and per-tract stat aggregation into a single TSV.
#
# Usage:
#   dwi_dsi_pipeline.sh [--skip-src] [--skip-rec] /path/to/sub-XXXX/dwi/
#
# Flags:
#   --skip-src    Skip SRC generation if .sz already exists
#   --skip-rec    Skip GQI reconstruction if .fib.gz already exists
#
# Expected inputs (outputs of dwi_fsl_preprocessing.sh):
#   eddy_unwarped_images.nii.gz
#   eddy_unwarped_images.eddy_rotated_bvecs
#   sub-*_dir-AP_dwi.bval
#   hifi_nodif_brain_mask.nii.gz        (optional; DSI Studio auto-mask used if absent)
#
# Outputs:
#   sub-XXXX_eddy_unwarped_images.sz
#   sub-XXXX_eddy_unwarped_images.fib.gz
#   tract_stats/
#     sub-XXXX_eddy_unwarped_images_autotrack.log
#     sub-XXXX_eddy_unwarped_images_autotrack_stats.tsv
#     sub-XXXX_eddy_unwarped_images_atk_work/
#
# === ADJUSTABLE PARAMETERS ===
#
# GQI reconstruction:
#   PARAM0            Diffusion sampling length ratio (default: 1.25 for in vivo).
#
# AutoTrack:
#   TOLERANCE         Bundle recognition tolerance in mm. Comma-separated values trigger
#                     progressive retry (default: "22,26,30" per CLI docs). Larger values
#                     accept more shape variation; increases false positives.
#   TRACK_VOXEL_RATIO Track-to-voxel ratio controlling streamline count (default: 2.0).
#                     Higher values improve yield but increase compute time.
#   CHECK_ENDING      Remove tracts terminating in high-anisotropy regions (default: 1).
#   YIELD_RATE        Early-termination threshold (default: 0.00001). Set to 0 to disable.
#   TIP_ITERATION     Topology-informed pruning iterations (default: 2).

PARAM0="1.25"

TOLERANCE="22,26,30"
TRACK_VOXEL_RATIO="2"
CHECK_ENDING="1"
YIELD_RATE="0.00001"
TIP_ITERATION="2"
THREAD_COUNT="${NSLOTS:-$(nproc)}"

# ==============================

SKIP_SRC=0
SKIP_REC=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-src) SKIP_SRC=1; shift ;;
    --skip-rec) SKIP_REC=1; shift ;;
    *) break ;;
  esac
done

module load dsi_studio

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 [--skip-src] [--skip-rec] /path/to/sub-XXXX/dwi/" >&2
  exit 1
fi

DWI_DIR="${1%/}"

if [ ! -d "$DWI_DIR" ]; then
  echo "ERROR: Directory not found: $DWI_DIR" >&2
  exit 1
fi

cd "$DWI_DIR"
echo "=== DSI PIPELINE: Working in $PWD ==="

# ---------- Check required inputs ----------

if [ ! -f "eddy_unwarped_images.nii.gz" ]; then
  echo "ERROR: eddy_unwarped_images.nii.gz not found in $PWD" >&2
  exit 1
fi

if [ ! -f "eddy_unwarped_images.eddy_rotated_bvecs" ]; then
  echo "ERROR: eddy_unwarped_images.eddy_rotated_bvecs not found in $PWD" >&2
  exit 1
fi

WM_MASK="hifi_nodif_brain_mask.nii.gz"
MASK_OPT=()
if [ -f "$WM_MASK" ]; then
  MASK_OPT=( "--mask=$WM_MASK" )
else
  echo "[WARN] $WM_MASK not found. DSI Studio will use its default mask."
fi

bval_file=$(ls sub-*_dir-AP_dwi.bval 2>/dev/null | head -n1 || true)
if [ -z "$bval_file" ]; then
  echo "ERROR: No sub-*_dir-AP_dwi.bval found in $PWD" >&2
  exit 1
fi

SUBJECT="${bval_file%%_dir-AP_dwi.bval}"
echo "[INFO] Subject: $SUBJECT"

SRC_OUT="${SUBJECT}_eddy_unwarped_images.sz"
FIB_OUT="${SUBJECT}_eddy_unwarped_images.fib.gz"

# ---------- Step 1: SRC ----------

if [[ "$SKIP_SRC" -eq 1 && -f "$SRC_OUT" ]]; then
  echo "[INFO] --skip-src set and SRC exists, skipping: $SRC_OUT"
else
  echo "[INFO] Generating SRC..."
  dsi_studio \
    --action=src \
    --source=eddy_unwarped_images.nii.gz \
    --bval="$bval_file" \
    --bvec=eddy_unwarped_images.eddy_rotated_bvecs \
    --output="$SRC_OUT"
  echo "[INFO] Created: $SRC_OUT"
fi

# ---------- Step 2: GQI reconstruction ----------

if [[ "$SKIP_REC" -eq 1 && -f "$FIB_OUT" ]]; then
  echo "[INFO] --skip-rec set and FIB exists, skipping: $FIB_OUT"
else
  echo "[INFO] Running GQI reconstruction..."
  dsi_studio \
    --action=rec \
    --source="$SRC_OUT" \
    --method=4 \
    --param0="${PARAM0}" \
    --other_output=fa,rd,md,ad,rdi \
    --thread_count="${THREAD_COUNT}" \
    "${MASK_OPT[@]}" \
    --output="$FIB_OUT"

  if [ ! -f "$FIB_OUT" ]; then
    ALT=$(ls "${SUBJECT}_eddy_unwarped_images"*.fib.gz 2>/dev/null | head -n1 || true)
    if [ -n "$ALT" ]; then
      echo "[INFO] Renaming $ALT -> $FIB_OUT"
      mv "$ALT" "$FIB_OUT"
    else
      echo "ERROR: Expected FIB output not found after reconstruction." >&2
      exit 1
    fi
  fi
  echo "[INFO] Created: $FIB_OUT"
fi

# ---------- Step 3: AutoTrack ----------

echo "[INFO] Parameters:"
echo "[INFO]   param0:            ${PARAM0}"
echo "[INFO]   tolerance:         ${TOLERANCE}"
echo "[INFO]   track_voxel_ratio: ${TRACK_VOXEL_RATIO}"
echo "[INFO]   check_ending:      ${CHECK_ENDING}"
echo "[INFO]   yield_rate:        ${YIELD_RATE}"
echo "[INFO]   tip_iteration:     ${TIP_ITERATION}"
echo "[INFO]   thread_count:      ${THREAD_COUNT}"

FIB_STEM="${FIB_OUT%.fib.gz}"
OUT_ROOT="${DWI_DIR}/tract_stats"
WORKDIR="${OUT_ROOT}/${FIB_STEM}_atk_work"
LOG_FILE="${OUT_ROOT}/${FIB_STEM}_autotrack.log"
FINAL_TSV="${OUT_ROOT}/${FIB_STEM}_autotrack_stats.tsv"

mkdir -p "${OUT_ROOT}"

if [[ -d "${WORKDIR}" ]]; then
  echo "ERROR: AutoTrack workdir already exists: ${WORKDIR}" >&2
  echo "       Move or remove it before re-running." >&2
  exit 1
fi

TMP_WORKDIR="$(mktemp -d "${OUT_ROOT}/${FIB_STEM}_atk_work_tmp_XXXXXX")"

echo "[INFO] Running AutoTrack..."
(
  dsi_studio \
    --action=atk \
    --source="${FIB_OUT}" \
    --trk_format=tt.gz \
    --tolerance="${TOLERANCE}" \
    --track_voxel_ratio="${TRACK_VOXEL_RATIO}" \
    --check_ending="${CHECK_ENDING}" \
    --yield_rate="${YIELD_RATE}" \
    --tip_iteration="${TIP_ITERATION}" \
    --thread_count="${THREAD_COUNT}" \
    --output="${TMP_WORKDIR}"
) 2>&1 | tee "${LOG_FILE}"

mv "${TMP_WORKDIR}" "${WORKDIR}"
echo "[INFO] AutoTrack complete."

# ---------- Step 4: Aggregate stats TSV ----------

echo "[INFO] Building final TSV..."

for f in "${WORKDIR}/${FIB_STEM}."*.tt.gz.stat.txt \
          "${WORKDIR}/${FIB_STEM}."*.no_result.txt \
          "${WORKDIR}/${FIB_STEM}."*.tt.gz; do
  [[ -e "$f" ]] || continue
  fname=$(basename "$f")
  tractname="${fname#${FIB_STEM}.}"
  tractname="${tractname%.tt.gz.stat.txt}"
  tractname="${tractname%.no_result.txt}"
  tractname="${tractname%.tt.gz}"
  mkdir -p "${WORKDIR}/${tractname}"
  mv "$f" "${WORKDIR}/${tractname}/"
done

mapfile -t TRACTS < <(
  cd "${WORKDIR}"
  LC_ALL=C printf '%s\n' */ | sed 's:/$::' | LC_ALL=C sort
)

if [[ ${#TRACTS[@]} -eq 0 ]]; then
  echo "ERROR: No tract subdirectories found in ${WORKDIR}" >&2
  exit 1
fi

FIRST_STAT=""
for t in "${TRACTS[@]}"; do
  shopt -s nullglob
  stats=( "${WORKDIR}/${t}"/*.stat.txt )
  shopt -u nullglob
  if [[ ${#stats[@]} -gt 0 ]]; then
    FIRST_STAT="${stats[0]}"
    break
  fi
done

if [[ -z "${FIRST_STAT}" ]]; then
  echo "ERROR: No *.stat.txt found in ${WORKDIR}" >&2
  echo "       If all tracts failed, only *.no_result.txt files may be present." >&2
  exit 1
fi

KEYS_FILE="$(mktemp)"
awk -F'\t' 'NF>=2 && $1!="" {print $1}' "${FIRST_STAT}" > "${KEYS_FILE}"

AGTMPDIR="$(mktemp -d)"
COL0="${AGTMPDIR}/col0.txt"
cp "${KEYS_FILE}" "${COL0}"

COLFILES=("${COL0}")
HEADER=("stat")

for t in "${TRACTS[@]}"; do
  HEADER+=("${t}")

  COL="${AGTMPDIR}/${t}.col.txt"
  shopt -s nullglob
  stats=( "${WORKDIR}/${t}"/*.stat.txt )
  shopt -u nullglob

  if [[ ${#stats[@]} -gt 0 ]]; then
    STATFILE="${stats[0]}"
    awk -F'\t' '
      NR==FNR { keys[++n]=$1; next }
      NF>=2 && $1!="" { val[$1]=$2 }
      END {
        for (i=1; i<=n; i++) {
          k=keys[i]
          if (k in val) print val[k]
          else print "NA"
        }
      }
    ' "${KEYS_FILE}" "${STATFILE}" > "${COL}"
  else
    awk '{print "NA"}' "${KEYS_FILE}" > "${COL}"
  fi

  COLFILES+=("${COL}")
done

{
  (IFS=$'\t'; echo "${HEADER[*]}")
  paste -d $'\t' "${COLFILES[@]}"
} > "${FINAL_TSV}"

rm -f "${KEYS_FILE}"
rm -rf "${AGTMPDIR}"

echo "[INFO] Wrote: ${FINAL_TSV}"
echo "[INFO] Workdir: ${WORKDIR}"
echo "=== DSI PIPELINE COMPLETE for ${SUBJECT} ==="
