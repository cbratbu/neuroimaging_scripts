#!/usr/bin/env bash
# build_conn_cfg_HC_single.sh
# Builds CONN .cfg files for healthy control subjects
# Optional flags: --useFLAIR, --useT2w to allow FLAIR/T2w as structural when T1w missing

set -euo pipefail

# ----- Adjustable fields -----
# Default TR (seconds) used if not found in JSON sidecar - update per cohort
DEFAULT_TR="2.0"
# Default slice order - update to match acquisition (e.g. "interleaved (Siemens)", "ascending", "descending")
DEFAULT_SLICEORDER="interleaved (Siemens)"

USE_FLAIR=0
USE_T2W=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --useFLAIR)
      USE_FLAIR=1
      shift
      ;;
    --useT2w)
      USE_T2W=1
      shift
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 [--useFLAIR] [--useT2w] /path/to/sub-XXXX" >&2
  exit 1
fi

SUBJ_DIR="${1%/}"
if [[ ! -d "$SUBJ_DIR" ]]; then
  echo "ERROR: subject dir not found: $SUBJ_DIR" >&2
  exit 1
fi

subj="$(basename "$SUBJ_DIR")"
if [[ "$subj" != sub-* ]]; then
  echo "ERROR: input must be a sub-* directory. Got: $SUBJ_DIR" >&2
  exit 1
fi

# Use 2.preprocessing if present; otherwise use subject root
BASE_DIR="$SUBJ_DIR"
if [[ -d "${SUBJ_DIR}/2.preprocessing" ]]; then
  BASE_DIR="${SUBJ_DIR}/2.preprocessing"
fi

CFG_DIR="${BASE_DIR}/cfg"
mkdir -p "$CFG_DIR"

shopt -s nullglob

pick_first_existing() {
  for f in "$@"; do
    [[ -f "$f" ]] && { echo "$f"; return 0; }
  done
  return 1
}

get_rt() {
  local func_path="$1"
  local json_path="${func_path%.nii.gz}"
  json_path="${json_path%.nii}"
  json_path="${json_path}.json"
  
  if [[ -f "$json_path" ]]; then
    rt=$(grep '"RepetitionTime"' "$json_path" | grep -o '[0-9.]\+' | head -1)
    if [[ -n "$rt" ]]; then
      echo "$rt"
      return 0
    fi
  fi
  return 1
}

sessions=("${BASE_DIR}"/ses-*/)

units=()
if (( ${#sessions[@]} > 0 )); then
  units=("${sessions[@]}")
else
  units=("${BASE_DIR}")
fi

echo "=== Building CONN cfgs (HC single-subject) ==="
echo "Subject dir : $SUBJ_DIR"
echo "Data base   : $BASE_DIR"
echo "CFG dir     : $CFG_DIR"
echo "Use FLAIR   : $USE_FLAIR"
echo "Use T2w     : $USE_T2W"
echo

wrote=0
skipped=0

for unit in "${units[@]}"; do
  unit="${unit%/}"

  ses=""
  cfg_id="$subj"
  if [[ "$(basename "$unit")" == ses-* ]]; then
    ses="$(basename "$unit")"
    cfg_id="${subj}_${ses}"
  fi

  # Build anatomical search list: T1w → FLAIR → T2w
  anat_candidates=()
  if [[ -n "$ses" ]]; then
    anat_candidates+=(
      "${unit}/anat/${subj}_${ses}_T1w.nii"
      "${unit}/anat/${subj}_${ses}_T1w.nii.gz"
    )
    if (( USE_FLAIR )); then
      anat_candidates+=(
        "${unit}/anat/${subj}_${ses}_FLAIR.nii"
        "${unit}/anat/${subj}_${ses}_FLAIR.nii.gz"
      )
    fi
    if (( USE_T2W )); then
      anat_candidates+=(
        "${unit}/anat/${subj}_${ses}_T2w.nii"
        "${unit}/anat/${subj}_${ses}_T2w.nii.gz"
      )
    fi
    
    func="$(pick_first_existing \
      "${unit}/func/${subj}_${ses}_task-rest_bold.nii" \
      "${unit}/func/${subj}_${ses}_task-rest_bold.nii.gz" \
    )" || func=""
  else
    anat_candidates+=(
      "${unit}/anat/${subj}_T1w.nii"
      "${unit}/anat/${subj}_T1w.nii.gz"
    )
    if (( USE_FLAIR )); then
      anat_candidates+=(
        "${unit}/anat/${subj}_FLAIR.nii"
        "${unit}/anat/${subj}_FLAIR.nii.gz"
      )
    fi
    if (( USE_T2W )); then
      anat_candidates+=(
        "${unit}/anat/${subj}_T2w.nii"
        "${unit}/anat/${subj}_T2w.nii.gz"
      )
    fi
    
    func="$(pick_first_existing \
      "${unit}/func/${subj}_task-rest_bold.nii" \
      "${unit}/func/${subj}_task-rest_bold.nii.gz" \
    )" || func=""
  fi

  anat="$(pick_first_existing "${anat_candidates[@]}")" || anat=""

  # Require structural only. Functional is optional.
  if [[ -z "$anat" ]]; then
    echo "[SKIP] ${cfg_id}: missing structural (T1w/FLAIR/T2w)"
    skipped=$((skipped+1))
    continue
  fi

  has_func=1
  rt="$DEFAULT_TR"
  
  if [[ -z "$func" ]]; then
    has_func=0
    echo "[NOFUNC] ${cfg_id}: missing task-rest BOLD; will write structural-only cfg"
  else
    rt=$(get_rt "$func") || rt="$DEFAULT_TR"
  fi

  cfg="${CFG_DIR}/${cfg_id}.cfg"

  if (( has_func )); then
    cat > "$cfg" <<EOF
#functionals
${func}

#structurals
${anat}

#steps
structural_center
structural_segment&normalize
functional_label_as_original
functional_realign
functional_center
functional_art
functional_label_as_subjectspace
functional_segment&normalize_direct
functional_label_as_mnispace
functional_regression
functional_label_as_denoised
functional_bandpass
functional_label_as_filtered
functional_smooth
functional_label_as_minimallysmoothed

#RT
${rt}

#sliceorder
${DEFAULT_SLICEORDER}

#fwhm
4

#reg_names
realignment
scrubbing
White Matter
CSF

#reg_dimensions
inf
inf
5
5

#reg_deriv
1
0
0
0

#bp_filter
0.008 0.09
EOF
  else
    cat > "$cfg" <<EOF
#functionals

#structurals
${anat}

#steps
structural_center
structural_segment&normalize

#RT
${rt}

#sliceorder
${DEFAULT_SLICEORDER}

#fwhm
4

#reg_names
realignment
scrubbing
White Matter
CSF

#reg_dimensions
inf
inf
5
5

#reg_deriv
1
0
0
0

#bp_filter
0.008 0.09
EOF
  fi

  echo "[OK] Wrote ${cfg} (RT=${rt})"
  wrote=$((wrote+1))
done

echo
echo "Done. Wrote: $wrote   Skipped: $skipped"