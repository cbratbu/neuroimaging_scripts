#!/usr/bin/env bash
# hcp_normative_ct_batch.sh
# Runs hcp_normative_ct_subject.sh in parallel across all valid HCP S1200
# subjects using GNU parallel, then aggregates results into normative CSV.
#
# Usage:
#   bash hcp_normative_ct_batch.sh [--jobs N] [--aggregate-only]
#
# Options:
#   --jobs N          Number of parallel jobs (default: 16)
#   --aggregate-only  Skip processing, just run aggregation on existing CSVs
#
# Load modules before running:
#   module load freesurfer/6.0 parallel

set -euo pipefail

# --- Adjustable fields -------------------------------------------------------
SCRIPT_DIR="$( cd "$(dirname "$0")" && pwd )"
HCP_BASE="/projectnb/connectomedb/HCP1200"
OUT_BASE="/projectnb/skiran/Computational_R01/3.Analysis/HCP1200_normative_CT"
PIPELINE_DIR="/projectnb/skiran/Computational_R01/2.Scripts/HO_CT_pipeline"
SUBJECT_SCRIPT="${SCRIPT_DIR}/hcp_normative_ct_subject.sh"
AGGREGATE_SCRIPT="${SCRIPT_DIR}/hcp_normative_ct_aggregate.py"
NORMS_OUT="${PIPELINE_DIR}/normative_thickness_ROI_stats.csv"
JOBLOG="${OUT_BASE}/parallel_joblog.txt"
N_JOBS=16

# --- Args --------------------------------------------------------------------
AGGREGATE_ONLY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --jobs)         N_JOBS="$2"; shift 2 ;;
    --aggregate-only) AGGREGATE_ONLY=1; shift ;;
    *) echo "WARNING: Unrecognized argument: $1"; shift ;;
  esac
done

mkdir -p "$OUT_BASE"

# --- Environment checks ------------------------------------------------------
[ -n "${FREESURFER_HOME:-}" ] || { echo "ERROR: FREESURFER_HOME not set. Load freesurfer/6.0"; exit 1; }
command -v parallel >/dev/null || { echo "ERROR: parallel not found. Run: module load parallel"; exit 1; }
command -v python3  >/dev/null || { echo "ERROR: python3 not found."; exit 1; }
[ -f "$SUBJECT_SCRIPT"   ] || { echo "ERROR: Missing: $SUBJECT_SCRIPT";   exit 1; }
[ -f "$AGGREGATE_SCRIPT" ] || { echo "ERROR: Missing: $AGGREGATE_SCRIPT"; exit 1; }

if [ "$AGGREGATE_ONLY" = "0" ]; then
  # --- Build subject list ----------------------------------------------------
  # Include only numeric subject IDs that have the expected FreeSurfer directory
  SUBJECT_LIST=$(mktemp)
  echo "[INFO] Scanning for valid subjects..."
  for ENTRY in "${HCP_BASE}"/*/; do
    ID="$(basename "$ENTRY")"
    [[ "$ID" =~ ^[0-9]+$ ]] || continue
    [ -d "${HCP_BASE}/${ID}/T1w/${ID}/surf" ] || continue
    echo "$ID"
  done > "$SUBJECT_LIST"

  N_TOTAL=$(wc -l < "$SUBJECT_LIST")
  echo "[INFO] Found ${N_TOTAL} valid subjects"
  echo "[INFO] Running with ${N_JOBS} parallel jobs"
  echo "[INFO] Job log: ${JOBLOG}"
  echo "[INFO] Resuming from joblog if interrupted (already-completed subjects skipped)"

  # GNU parallel:
  #   --joblog: tracks completed jobs; --resume-failed reruns only failed/unstarted
  #   --env FREESURFER_HOME: passes env var to each job
  #   --delay 0.2: stagger starts slightly to avoid simultaneous file I/O at startup
  #   -j N_JOBS: parallel workers
  parallel \
    --joblog "$JOBLOG" \
    --resume-failed \
    --env FREESURFER_HOME \
    --env PATH \
    --delay 0.2 \
    -j "$N_JOBS" \
    bash "$SUBJECT_SCRIPT" {} "$HCP_BASE" "$OUT_BASE" \
    :::: "$SUBJECT_LIST"

  rm -f "$SUBJECT_LIST"

  # Count completed
  N_DONE=$(find "$OUT_BASE" -name "*_HO48_thickness.csv" -not -path "*/stub/*" | wc -l)
  echo "[INFO] Completed: ${N_DONE}/${N_TOTAL} subjects"
fi

# --- Aggregation -------------------------------------------------------------
echo "[INFO] Running aggregation..."
python3 "$AGGREGATE_SCRIPT" "$OUT_BASE" "$NORMS_OUT"
echo "[DONE] Normative CSV: $NORMS_OUT"
