#!/usr/bin/env bash
# batch_run_conn_on_cfg.sh - with resumability at subject level
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/parent_dir" >&2
  exit 1
fi

PARENT_DIR="${1%/}"
if [[ ! -d "$PARENT_DIR" ]]; then
  echo "ERROR: parent dir not found: $PARENT_DIR" >&2
  exit 1
fi

SCRIPT_DIR="/projectnb/openneuro-aphasia/ML_aphasia_scripts/CONN_Preprocessing"
RUNNER="${SCRIPT_DIR}/run_conn_on_cfg.sh"

if [[ ! -x "$RUNNER" ]]; then
  echo "ERROR: runner not found or not executable: $RUNNER" >&2
  exit 1
fi

mapfile -t SUBJECTS < <(find "$PARENT_DIR" -maxdepth 1 -type d -name "sub-*" | sort)
TOTAL=${#SUBJECTS[@]}

echo "=== Batch run CONN preprocessing ==="
echo "Parent dir : ${PARENT_DIR}"
echo "Runner     : ${RUNNER}"
echo "Found      : ${TOTAL} subjects"
echo

ok_subjects=0
skip_subjects=0
fail_subjects=0

failed_list=()
skipped_list=()

i=0
for SUBJ_DIR in "${SUBJECTS[@]}"; do
  i=$((i+1))
  subj="$(basename "$SUBJ_DIR")"
  
  # Check if subject already processed by looking at all cfgs
  BASE_DIR="$SUBJ_DIR"
  if [[ -d "${SUBJ_DIR}/2.preprocessing" ]]; then
    BASE_DIR="${SUBJ_DIR}/2.preprocessing"
  fi
  
  CFG_DIR="${BASE_DIR}/cfg"
  
  if [[ ! -d "$CFG_DIR" ]]; then
    echo "=== [${i}/${TOTAL}] ${subj} ==="
    echo "[SKIP] ${subj}: no cfg directory"
    skip_subjects=$((skip_subjects+1))
    skipped_list+=("$subj")
    echo
    continue
  fi
  
  mapfile -t CFGS < <(find "$CFG_DIR" -maxdepth 1 -type f -name "*.cfg" 2>/dev/null | sort)
  
  if (( ${#CFGS[@]} == 0 )); then
    echo "=== [${i}/${TOTAL}] ${subj} ==="
    echo "[SKIP] ${subj}: no cfgs found"
    skip_subjects=$((skip_subjects+1))
    skipped_list+=("$subj")
    echo
    continue
  fi
  
  # Check if all cfgs are already processed
  all_done=true
  for CFG in "${CFGS[@]}"; do
    func_path=$(awk '/^#functionals$/ {getline; print; exit}' "$CFG")
    
    if [[ -n "$func_path" && -f "$func_path" ]]; then
      # Has functional data - check for final smoothed output
      func_dir=$(dirname "$func_path")
      func_base=$(basename "$func_path")
      func_base="${func_base%.nii.gz}"
      func_base="${func_base%.nii}"
      
      final_func="${func_dir}/sbdmwr${func_base}.nii"
      
      if [[ ! -f "$final_func" ]]; then
        all_done=false
        break
      fi
    else
      # Structural-only - check for final structural output
      struct_path=$(awk '/^#structurals$/ {getline; print; exit}' "$CFG")
      
      if [[ -n "$struct_path" && -f "$struct_path" ]]; then
        struct_dir=$(dirname "$struct_path")
        struct_base=$(basename "$struct_path")
        struct_base="${struct_base%.nii.gz}"
        struct_base="${struct_base%.nii}"
        
        final_struct="${struct_dir}/wc0m${struct_base}.nii"
        
        if [[ ! -f "$final_struct" ]]; then
          all_done=false
          break
        fi
      else
        all_done=false
        break
      fi
    fi
  done
  
  if $all_done; then
    echo "=== [${i}/${TOTAL}] ${subj} ==="
    echo "[SKIP] ${subj}: all cfgs already processed"
    skip_subjects=$((skip_subjects+1))
    skipped_list+=("$subj")
    echo
    continue
  fi

  echo "=== [${i}/${TOTAL}] ${subj} ==="

  out="$("$RUNNER" "$SUBJ_DIR" 2>&1)" || {
    # If no cfgs, treat as SKIP instead of FAIL
    if printf "%s\n" "$out" | grep -qE 'ERROR: (cfg directory not found|no \.cfg files found)'; then
      echo "$out"
      echo "[SKIP] ${subj}: no cfgs found"
      skip_subjects=$((skip_subjects+1))
      skipped_list+=("$subj")
      echo
      continue
    fi

    echo "$out"
    echo "[FAIL] ${subj}"
    fail_subjects=$((fail_subjects+1))
    failed_list+=("$subj")
    echo
    continue
  }

  echo "$out"
  echo "[OK] ${subj}"
  ok_subjects=$((ok_subjects+1))
  echo
done

echo "=== Batch CONN run summary ==="
echo "Subjects total  : $TOTAL"
echo "Subjects OK     : $ok_subjects"
echo "Subjects SKIP   : $skip_subjects"
echo "Subjects FAIL   : $fail_subjects"

if (( skip_subjects > 0 )); then
  echo "Skipped subjects: ${skipped_list[*]}"
fi
if (( fail_subjects > 0 )); then
  echo "Failed subjects : ${failed_list[*]}"
fi

(( fail_subjects == 0 )) || exit 1
