#!/usr/bin/env bash
# run_conn_on_cfg.sh
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /path/to/sub-XXXX" >&2
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

module purge
module load matlab/2023a spm/12.7771 conn/

if [[ ! -d "$CFG_DIR" ]]; then
  echo "ERROR: cfg directory not found: $CFG_DIR" >&2
  exit 1
fi

mapfile -t CFGS < <(find "$CFG_DIR" -maxdepth 1 -type f -name "*.cfg" | sort)
TOTAL=${#CFGS[@]}

if (( TOTAL == 0 )); then
  echo "ERROR: no .cfg files found in: $CFG_DIR" >&2
  exit 1
fi

echo "=== Running CONN preprocessing for ${subj} ==="
echo "Subject dir : $SUBJ_DIR"
echo "Data base   : $BASE_DIR"
echo "CFG dir     : $CFG_DIR"
echo "Found cfgs  : $TOTAL"
printf '  %s\n' "${CFGS[@]}"
echo

ok=0
fail=0
i=0

for CFG in "${CFGS[@]}"; do
  i=$((i+1))
  id="$(basename "${CFG%.*}")"
  echo "=== [${i}/${TOTAL}] Running ${id} ==="
  echo "CFG : ${CFG}"

  if matlab -nodisplay -nosplash -r "try; addpath('/share/pkg.8/conn/22v2407/install/conn'); addpath('/share/pkg.7/spm/12.7771/install/spm12'); conn_module('preprocessing','${CFG}'); exit(0); catch ME; disp(getReport(ME)); exit(1); end"; then
    echo "[OK] ${id}"
    ok=$((ok+1))
  else
    echo "[FAIL] ${id}"
    fail=$((fail+1))
  fi
  echo
done

echo "=== CONN summary: ${subj} ==="
echo "Total cfgs: $TOTAL"
echo "OK       : $ok"
echo "Fail     : $fail"
(( fail == 0 )) || exit 1
