#!/usr/bin/env bash
# hcp_normative_ct_subject.sh
# Processes one HCP S1200 subject: transfers HO thr25 atlas labels to native
# surface via mri_label2label, assembles annotation, runs mris_anatomical_stats,
# writes per-subject thickness and volumes CSVs.
# Called by hcp_normative_ct_batch.sh — can also be run standalone.
#
# Usage:
#   bash hcp_normative_ct_subject.sh <subject_id> <hcp_base_dir> <output_base_dir>
#
# Example:
#   bash hcp_normative_ct_subject.sh 100206 \
#     /projectnb/connectomedb/HCP1200 \
#     /projectnb/skiran/Computational_R01/3.Analysis/HCP1200_normative_CT
#
# Load modules before running:
#   module load freesurfer/6.0

set -euo pipefail

# --- Args --------------------------------------------------------------------
if [ $# -lt 3 ]; then
  echo "Usage: $0 <subject_id> <hcp_base_dir> <output_base_dir>"
  exit 1
fi

SUB="$1"
HCP_BASE="$2"
OUT_BASE="$3"

# --- Globals -----------------------------------------------------------------
SCRIPT_DIR="$( cd "$(dirname "$0")" && pwd )"
PIPELINE_DIR="/projectnb/skiran/Computational_R01/2.Scripts/HO_CT_pipeline"
ATLAS_DIR="${PIPELINE_DIR}/atlas"
CTAB="${SCRIPT_DIR}/HO48_1based.lut"
LABELS_TSV="${SCRIPT_DIR}/HO48_labels_1based.tsv"
FSAVG_LABEL_DIR="${ATLAS_DIR}/fsaverage_labels/thr25"
FSAVERAGE_SRC="${ATLAS_DIR}/subjects_tmp/fsaverage"
THR="thr25"

# HCP FreeSurfer outputs live at HCP_BASE/<subid>/T1w/<subid>/
HCP_SUB_DIR="${HCP_BASE}/${SUB}/T1w/${SUB}"

# --- Validate subject --------------------------------------------------------
if [ ! -d "${HCP_SUB_DIR}/surf" ]; then
  echo "[SKIP] ${SUB}: FreeSurfer dir not found at ${HCP_SUB_DIR}"
  exit 0
fi
for f in \
  "${HCP_SUB_DIR}/surf/lh.sphere.reg" \
  "${HCP_SUB_DIR}/surf/rh.sphere.reg" \
  "${HCP_SUB_DIR}/surf/lh.thickness" \
  "${HCP_SUB_DIR}/surf/rh.thickness" \
  "${HCP_SUB_DIR}/surf/lh.white" \
  "${HCP_SUB_DIR}/surf/rh.white" \
  "${HCP_SUB_DIR}/surf/lh.pial" \
  "${HCP_SUB_DIR}/surf/rh.pial" \
  "${HCP_SUB_DIR}/surf/lh.sulc" \
  "${HCP_SUB_DIR}/surf/rh.sulc" \
  "${HCP_SUB_DIR}/label/lh.cortex.label" \
  "${HCP_SUB_DIR}/label/rh.cortex.label" \
  "${HCP_SUB_DIR}/mri/ribbon.mgz" \
  "${HCP_SUB_DIR}/stats/aseg.stats"; do
  if [ ! -f "$f" ]; then
    echo "[SKIP] ${SUB}: Missing required file: $f"
    exit 0
  fi
done

# --- Output paths ------------------------------------------------------------
SUB_OUT="${OUT_BASE}/${SUB}"
FINAL_CSV="${SUB_OUT}/${SUB}_HO48_thickness.csv"
VOLCSV="${SUB_OUT}/${SUB}_subject_volumes.csv"

# Skip if already complete
if [ -f "$FINAL_CSV" ] && [ -f "$VOLCSV" ]; then
  echo "[SKIP] ${SUB}: Already complete"
  exit 0
fi

mkdir -p "$SUB_OUT"

# --- Environment checks ------------------------------------------------------
[ -n "${FREESURFER_HOME:-}" ] || { echo "ERROR: FREESURFER_HOME not set."; exit 1; }
for cmd in mri_label2label mris_label2annot mris_anatomical_stats python3; do
  command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not found."; exit 1; }
done

# --- Build subject stub ------------------------------------------------------
# mri_label2label and mris_label2annot require SUBJECTS_DIR to contain both
# fsaverage and the subject. We create a stub with a real writable label/ dir
# and symlinks to the read-only HCP directories for everything else.
# mris_label2annot writes the annot to stub/label/ — nothing touches HCP dir.
STUB="${SUB_OUT}/stub"
mkdir -p "${STUB}/label"
[ -e "${STUB}/surf"     ] || ln -s "${HCP_SUB_DIR}/surf"     "${STUB}/surf"
[ -e "${STUB}/mri"      ] || ln -s "${HCP_SUB_DIR}/mri"      "${STUB}/mri"
[ -e "${STUB}/stats"    ] || ln -s "${HCP_SUB_DIR}/stats"    "${STUB}/stats"

# Symlink all existing HCP label files into stub/label/ so mris_anatomical_stats
# can find cortex.label and any other label files it needs. New files written
# by mris_label2annot land in the same real directory without conflict.
for LFILE in "${HCP_SUB_DIR}/label"/*; do
  LNAME="$(basename "$LFILE")"
  [ -e "${STUB}/label/${LNAME}" ] || ln -s "$LFILE" "${STUB}/label/${LNAME}"
done

FS_WORK="${SUB_OUT}/subjects_work"
mkdir -p "$FS_WORK"
[ -e "${FS_WORK}/fsaverage" ] || ln -s "$(cd "$FSAVERAGE_SRC" && pwd)" "${FS_WORK}/fsaverage"
[ -e "${FS_WORK}/${SUB}"    ] || ln -s "$(cd "$STUB" && pwd)"          "${FS_WORK}/${SUB}"
export SUBJECTS_DIR="$FS_WORK"

trap 'rm -f "${FS_WORK}/fsaverage" "${FS_WORK}/${SUB}"' EXIT

# --- Step 1: Transfer HO labels from fsaverage to subject surface ------------
LABEL_DIR="${SUB_OUT}/labels"
mkdir -p "$LABEL_DIR"

for HEMI in lh rh; do
  ALL_PRESENT=1
  for IDX in $(seq 1 48); do
    [ -f "${LABEL_DIR}/${HEMI}.ho_${IDX}.label" ] || { ALL_PRESENT=0; break; }
  done
  [ "$ALL_PRESENT" = "1" ] && continue

  echo "[RUN] ${SUB}: mri_label2label ${HEMI}"
  for IDX in $(seq 1 48); do
    OUT_LABEL="${LABEL_DIR}/${HEMI}.ho_${IDX}.label"
    [ -f "$OUT_LABEL" ] && continue
    mri_label2label \
      --srcsubject fsaverage \
      --srclabel   "${FSAVG_LABEL_DIR}/${HEMI}/${HEMI}.HO_${IDX}.label" \
      --trgsubject "$SUB" \
      --trglabel   "$OUT_LABEL" \
      --hemi       "$HEMI" \
      --regmethod  surface
  done
done

# --- Step 2: Assemble subject-space annotation -------------------------------
for HEMI in lh rh; do
  ANNOT_OUT="${SUB_OUT}/${HEMI}.ho_${THR}.annot"
  [ -f "$ANNOT_OUT" ] && continue

  echo "[RUN] ${SUB}: mris_label2annot ${HEMI}"
  LABEL_ARGS=""
  for IDX in $(seq 1 48); do
    LFILE="${LABEL_DIR}/${HEMI}.ho_${IDX}.label"
    [ -f "$LFILE" ] && LABEL_ARGS="$LABEL_ARGS --l $LFILE"
  done

  # shellcheck disable=SC2086
  mris_label2annot \
    --sd   "$SUBJECTS_DIR" \
    --s    "$SUB" \
    --ctab "$CTAB" \
    $LABEL_ARGS \
    --h    "$HEMI" \
    --a    "ho_${THR}"

  # mris_label2annot writes to stub/label/; copy to output dir
  SRC="${STUB}/label/${HEMI}.ho_${THR}.annot"
  [ -f "$SRC" ] || { echo "ERROR: ${SUB}: Annot not found: $SRC"; exit 1; }
  cp "$SRC" "$ANNOT_OUT"
done

# --- Step 3: mris_anatomical_stats per hemisphere ----------------------------
for HEMI in lh rh; do
  STATS_OUT="${SUB_OUT}/${HEMI}.ho_${THR}.stats"
  [ -f "$STATS_OUT" ] && continue

  echo "[RUN] ${SUB}: mris_anatomical_stats ${HEMI}"
  mris_anatomical_stats \
    -mgz \
    -cortex "${HCP_SUB_DIR}/label/${HEMI}.cortex.label" \
    -f      "$STATS_OUT" \
    -b \
    -a      "${SUB_OUT}/${HEMI}.ho_${THR}.annot" \
    -c      "$CTAB" \
    "$SUB" "$HEMI"
done

# --- Step 4: Parse stats to CSV ----------------------------------------------
LH_STATS="${SUB_OUT}/lh.ho_${THR}.stats"
RH_STATS="${SUB_OUT}/rh.ho_${THR}.stats"
LABELS_TSV_P="$LABELS_TSV"
CTAB_P="$CTAB"
FINAL_CSV_P="$FINAL_CSV"
SUB_P="$SUB"

python3 - <<PY
import csv, math

lh_stats  = "$LH_STATS"
rh_stats  = "$RH_STATS"
tsv_path  = "$LABELS_TSV_P"
ctab_path = "$CTAB_P"
out_path  = "$FINAL_CSV_P"
sub       = "$SUB_P"

idx2ctab = {}
idx2display = {}
with open(ctab_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        p = line.split()
        if len(p) >= 2: idx2ctab[int(p[0])] = p[1]
with open(tsv_path) as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) >= 2: idx2display[int(p[0])] = p[1]

def parse_stats(path):
    data = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            p = line.split()
            if len(p) < 6: continue
            data[p[0]] = {'nvox': int(p[1]), 'vol_mm3': float(p[3]),
                          'mean': float(p[4]), 'sd': float(p[5])}
    return data

def fmt(x): return "NaN" if x != x else f"{x:.6f}"
def nanf(x):
    try: return float(x)
    except: return float('nan')

lh = parse_stats(lh_stats)
rh = parse_stats(rh_stats)

fields = ['SubjectID','Index','ROIName',
    'LH_mean','LH_sd','LH_nvox','LH_vol_mm3',
    'RH_mean','RH_sd','RH_nvox','RH_vol_mm3',
    'Mean_CT','LI']

out_rows = []
for idx in range(1, 49):
    name    = idx2ctab.get(idx, f"ROI_{idx}")
    display = idx2display.get(idx, name)
    lh_d = lh.get(name, {})
    rh_d = rh.get(name, {})
    lm = nanf(lh_d.get('mean', 'nan') if lh_d else 'nan')
    rm = nanf(rh_d.get('mean', 'nan') if rh_d else 'nan')
    if lm==lm and rm==rm: mct=(lm+rm)/2; li=(lm-rm)/(lm+rm) if (lm+rm)!=0 else float('nan')
    elif lm==lm:           mct=lm;        li=float('nan')
    elif rm==rm:           mct=rm;        li=float('nan')
    else:                  mct=float('nan'); li=float('nan')
    out_rows.append({
        'SubjectID': sub, 'Index': idx, 'ROIName': display,
        'LH_mean': fmt(lm), 'LH_sd': fmt(nanf(lh_d.get('sd','nan') if lh_d else 'nan')),
        'LH_nvox': lh_d.get('nvox',0) if lh_d else 0,
        'LH_vol_mm3': fmt(nanf(lh_d.get('vol_mm3','nan') if lh_d else 'nan')),
        'RH_mean': fmt(rm), 'RH_sd': fmt(nanf(rh_d.get('sd','nan') if rh_d else 'nan')),
        'RH_nvox': rh_d.get('nvox',0) if rh_d else 0,
        'RH_vol_mm3': fmt(nanf(rh_d.get('vol_mm3','nan') if rh_d else 'nan')),
        'Mean_CT': fmt(mct), 'LI': fmt(li),
    })

with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(out_rows)

n_lh = sum(1 for r in out_rows if r['LH_mean'] != 'NaN')
n_rh = sum(1 for r in out_rows if r['RH_mean'] != 'NaN')
print(f"[OK] {sub}: LH {n_lh}/48  RH {n_rh}/48 -> {out_path}")
PY

# --- Step 5: Subject volumes from aseg.stats ---------------------------------
ASEG="${HCP_SUB_DIR}/stats/aseg.stats"
VOLCSV_P="$VOLCSV"
SUB_P2="$SUB"

python3 - <<PY
import re, csv
aseg_path = "$ASEG"
out_path  = "$VOLCSV_P"
sub       = "$SUB_P2"

MEASURE_PAT = re.compile(r'^\s*#\s*Measure\s+([^,]+)\s*,\s*([^,]+)\s*,[^,]*,\s*([0-9eE.+-]+)\s*,')
VENTRICLE_STRUCTS = {'Left-Lateral-Ventricle','Right-Lateral-Ventricle','3rd-Ventricle',
                     '4th-Ventricle','Left-choroid-plexus','Right-choroid-plexus'}
vals = {}
vent_vol = 0.0
vent_found = set()
with open(aseg_path, 'r', errors='ignore') as f:
    for line in f:
        m = MEASURE_PAT.match(line)
        if m:
            vals[m.group(1).strip()] = float(m.group(3))
            vals[m.group(2).strip()] = float(m.group(3))
            continue
        if line.startswith('#'): continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                if parts[4] in VENTRICLE_STRUCTS:
                    vent_vol += float(parts[3])
                    vent_found.add(parts[4])
            except (ValueError, IndexError): pass

def get(*keys):
    for k in keys:
        if k in vals: return vals[k]
    return None
def fmt(x):
    try: return f"{float(x):.6f}"
    except: return "NaN"

eTIV = get('EstimatedTotalIntraCranialVol','eTIV')
GM   = get('TotalGrayVol','TotalGray')
WM   = get('CerebralWhiteMatterVol','CerebralWhiteMatter','CorticalWhiteMatterVol','CorticalWhiteMatter')
BSV  = get('BrainSegVol','BrainSeg')
CSF  = (BSV-(GM+WM)) if all(x is not None for x in [BSV,GM,WM]) else None
VC   = vent_vol if vent_found else None

with open(out_path,'w',newline='') as f:
    w = csv.writer(f)
    w.writerow(['SubjectID','Measure','Value_mm3'])
    for name,val in [
        ('EstimatedTotalIntraCranialVol',eTIV),
        ('TotalGrayVol',GM),
        ('CerebralWhiteMatterVol',WM),
        ('CSF_total',CSF),
        ('CSF_VentricleChoroidVol',VC)]:
        w.writerow([sub, name, fmt(val)])
print(f"[OK] {sub}: volumes -> {out_path}")
PY

echo "[DONE] ${SUB}"
