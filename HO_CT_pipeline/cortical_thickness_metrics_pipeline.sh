#!/usr/bin/env bash
# cortical_thickness_metrics_pipeline.sh
# Per-subject cortical thickness extraction in Harvard-Oxford ROIs using
# FreeSurfer surface-based registration (sphere.reg). HO atlas labels are
# transferred from fsaverage to the subject's native surface via
# mri_label2label, then mris_anatomical_stats extracts mean thickness per ROI.
# No volumetric registration to MNI space is performed.
#
# Usage:
#   bash cortical_thickness_metrics_pipeline.sh /path/to/<freesurfer_subject_dir> \
#     [out=/path/to/output_dir] [thr=thr25|thr0]
#
# Adjustable fields:
#   GLOBAL_DIR  — pipeline reference directory
#   DEFAULT_THR — atlas threshold variant (thr25 recommended)
#
# Load modules before running:
#   module load freesurfer/6.0

set -euo pipefail

DEFAULT_THR="thr25"

# --- Args --------------------------------------------------------------------
if [ $# -lt 1 ]; then
  echo "Usage: $0 /path/to/<freesurfer_subject_dir> [out=/path/to/output_dir] [thr=thr25|thr0]"
  exit 1
fi

SUBDIR="$1"
OUTDIR=""
THR="$DEFAULT_THR"

for arg in "${@:2}"; do
  case "$arg" in
    out=*) OUTDIR="${arg#out=}" ;;
    thr=*) THR="${arg#thr=}"   ;;
    *)     echo "WARNING: Unrecognized argument ignored: $arg" ;;
  esac
done

[ -d "$SUBDIR/mri" ] || { echo "ERROR: Not a FreeSurfer subject dir: $SUBDIR"; exit 1; }
[[ "$THR" == "thr0" || "$THR" == "thr25" ]] || {
  echo "ERROR: thr must be thr0 or thr25, got: $THR"; exit 1; }

[ -z "$OUTDIR" ] && OUTDIR="$SUBDIR"
BASE_OUTDIR="$OUTDIR"

# --- Globals -----------------------------------------------------------------
GLOBAL_DIR="/projectnb/skiran/Computational_R01/2.Scripts/HO_CT_pipeline"
ATLAS_DIR="${GLOBAL_DIR}/atlas"
CTAB="${GLOBAL_DIR}/HO48_1based.lut"
LABELS_TSV="${GLOBAL_DIR}/HO48_labels_1based.tsv"
NORMS_CSV="${GLOBAL_DIR}/normative_thickness_ROI_stats.csv"
VOLNORMS_CSV="${GLOBAL_DIR}/normative_subject_volume_means.csv"

FSAVG_LABEL_DIR="${ATLAS_DIR}/fsaverage_labels/${THR}"
FSAVERAGE_SRC="${ATLAS_DIR}/subjects_tmp/fsaverage"

for f in "$CTAB" "$LABELS_TSV"; do
  [ -f "$f" ] || { echo "ERROR: Missing file: $f"; exit 1; }
done
for hemi in lh rh; do
  [ -f "${ATLAS_DIR}/${hemi}.ho_${THR}.annot" ] || {
    echo "ERROR: Missing atlas annot: ${ATLAS_DIR}/${hemi}.ho_${THR}.annot"
    echo "       Run ho_ct_setup_fsaverage.sh first."; exit 1; }
  for idx in $(seq 1 48); do
    [ -f "${FSAVG_LABEL_DIR}/${hemi}/${hemi}.HO_${idx}.label" ] || {
      echo "ERROR: Missing setup label: ${FSAVG_LABEL_DIR}/${hemi}/${hemi}.HO_${idx}.label"
      echo "       Run ho_ct_setup_fsaverage.sh first."; exit 1; }
  done
done
[ -d "${FSAVERAGE_SRC}/surf" ] || {
  echo "ERROR: Writable fsaverage not found at ${FSAVERAGE_SRC}"
  echo "       Run ho_ct_setup_fsaverage.sh first."; exit 1; }

# --- Environment checks ------------------------------------------------------
echo "[INFO] FREESURFER_HOME=${FREESURFER_HOME:-<unset>}"
[ -n "${FREESURFER_HOME:-}" ] || { echo "ERROR: FREESURFER_HOME not set."; exit 1; }
for cmd in mri_label2label mris_label2annot mris_anatomical_stats python3; do
  command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not found on PATH."; exit 1; }
done

# --- Subject identity --------------------------------------------------------
SUB="$(basename "$SUBDIR")"
STATS_DIR="$SUBDIR/stats"
SUB_LABEL_DIR="$SUBDIR/label"

SUB_ID=$(echo "$SUBDIR" | grep -oE 'sub-[^/]+' | head -n 1 || true)
SES_ID=$(echo "$SUBDIR" | grep -oE 'ses-[^/]+' | head -n 1 || true)
[ -z "$SUB_ID" ] && SUB_ID="$SUB"
[ -n "$SES_ID" ] && OUT_PREFIX="${SUB_ID}_${SES_ID}" || OUT_PREFIX="${SUB_ID}"

# All outputs written into a named subfolder. Delete this folder to rerun cleanly.
OUTDIR="${BASE_OUTDIR}/${OUT_PREFIX}_cortical_thickness_pipeline"
mkdir -p "$OUTDIR"

WORKDIR="$OUTDIR/work"
mkdir -p "$WORKDIR"

FINAL_CSV="$OUTDIR/${OUT_PREFIX}_HO48_thickness.csv"
VOLCSV="$OUTDIR/${OUT_PREFIX}_subject_volumes.csv"
PVSN_CSV="$OUTDIR/${OUT_PREFIX}_patient_vs_norms.csv"
VOLPVSN_CSV="$OUTDIR/${OUT_PREFIX}_volumes_vs_norms.csv"

echo "[INFO] Subject:    $SUB"
echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] Atlas thr:  $THR"

# --- SUBJECTS_DIR setup ------------------------------------------------------
# FreeSurfer tools require both fsaverage and the subject to be present under
# a common SUBJECTS_DIR. We create a lightweight per-run SUBJECTS_DIR inside
# WORKDIR containing only symlinks — no data is copied.
# mris_label2annot writes the subject annot to SUBDIR/label/ (via the symlink),
# which is the correct location and stays with the subject's recon-all outputs.
# Nothing is written to the global pipeline directory.
FS_WORK="${WORKDIR}/subjects"
mkdir -p "$FS_WORK"
[ -e "${FS_WORK}/fsaverage" ] || ln -s "$(cd "$FSAVERAGE_SRC" && pwd)" "${FS_WORK}/fsaverage"
[ -e "${FS_WORK}/${SUB}"    ] || ln -s "$(cd "$SUBDIR" && pwd)"        "${FS_WORK}/${SUB}"
export SUBJECTS_DIR="$FS_WORK"

# Cleanup symlinks on exit (normal or error) so WORKDIR stays portable
trap 'rm -f "${FS_WORK}/fsaverage" "${FS_WORK}/${SUB}"' EXIT

# --- Step 1: Transfer HO labels from fsaverage to subject native surface -----
# mri_label2label --regmethod surface uses the subject's surf/lh.sphere.reg
# (produced by recon-all) to transfer labels via sulcal topology alignment.
# Transferred labels are stored in WORKDIR/labels/ to keep them separate from
# the subject's recon-all label outputs.
TRANSFERRED_LABEL_DIR="$WORKDIR/labels"
mkdir -p "$TRANSFERRED_LABEL_DIR"

for HEMI in lh rh; do
  ALL_PRESENT=1
  for IDX in $(seq 1 48); do
    [ -f "${TRANSFERRED_LABEL_DIR}/${HEMI}.ho_${IDX}.label" ] || { ALL_PRESENT=0; break; }
  done
  if [ "$ALL_PRESENT" = "1" ]; then
    echo "[SKIP] Labels already transferred: $HEMI"
    continue
  fi

  echo "[RUN] mri_label2label: fsaverage -> $SUB ($HEMI, 48 labels)"
  for IDX in $(seq 1 48); do
    OUT_LABEL="${TRANSFERRED_LABEL_DIR}/${HEMI}.ho_${IDX}.label"
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

# --- Step 2: Assemble subject-space HO annotation ----------------------------
# mris_label2annot writes to $SUBJECTS_DIR/$SUB/label/ (i.e. SUBDIR/label/).
# The annot is also copied to WORKDIR for use by mris_anatomical_stats and QC.
for HEMI in lh rh; do
  ANNOT_SUBDIR="${SUB_LABEL_DIR}/${HEMI}.ho_${THR}.annot"
  ANNOT_WORK="${WORKDIR}/${HEMI}.ho_${THR}.annot"

  if [ -f "$ANNOT_WORK" ]; then
    echo "[SKIP] Annotation exists: $ANNOT_WORK"
    continue
  fi

  echo "[RUN] mris_label2annot: $SUB $HEMI"
  # Remove any existing subject-space annot so mris_label2annot can write fresh.
  # This handles reruns after the output folder has been deleted.
  rm -f "${SUB_LABEL_DIR}/${HEMI}.ho_${THR}.annot"
  LABEL_ARGS=""
  for IDX in $(seq 1 48); do
    LFILE="${TRANSFERRED_LABEL_DIR}/${HEMI}.ho_${IDX}.label"
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

  # mris_label2annot writes to SUBDIR/label/ via the symlink
  [ -f "$ANNOT_SUBDIR" ] || { echo "ERROR: Annot not found: $ANNOT_SUBDIR"; exit 1; }
  cp "$ANNOT_SUBDIR" "$ANNOT_WORK"
done

# --- Step 3: Extract per-ROI thickness with mris_anatomical_stats ------------
# -cortex restricts sampling to cortical vertices (excludes medial wall).
# Stats are written to WORKDIR and parsed in Step 4.
for HEMI in lh rh; do
  STATS_OUT="${WORKDIR}/${HEMI}.ho_${THR}.stats"
  if [ -f "$STATS_OUT" ]; then
    echo "[SKIP] Stats exist: $STATS_OUT"
    continue
  fi
  echo "[RUN] mris_anatomical_stats: $SUB $HEMI"
  mris_anatomical_stats \
    -mgz \
    -cortex "${SUB_LABEL_DIR}/${HEMI}.cortex.label" \
    -f      "$STATS_OUT" \
    -b \
    -a      "${WORKDIR}/${HEMI}.ho_${THR}.annot" \
    -c      "$CTAB" \
    "$SUB" "$HEMI"
done

# --- Step 4: Parse stats and write thickness CSV -----------------------------
LH_STATS="$WORKDIR/lh.ho_${THR}.stats"
RH_STATS="$WORKDIR/rh.ho_${THR}.stats"
LABELS_TSV_PASS="$LABELS_TSV"
FINAL_CSV_PASS="$FINAL_CSV"

python3 - <<PY
import csv, math, os

lh_stats = "$LH_STATS"
rh_stats = "$RH_STATS"
tsv_path = "$LABELS_TSV_PASS"
out_path = "$FINAL_CSV_PASS"

# Build index->name from ctab (underscore names matching the stats file),
# and a separate display name map from the TSV (human-readable, for ROIName column).
idx2ctab = {}  # e.g. "Frontal_Pole"
idx2display = {}  # e.g. "Frontal Pole"
ctab_path = "$CTAB"
with open(ctab_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            idx2ctab[int(parts[0])] = parts[1]
with open(tsv_path) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            idx2display[int(parts[0])] = parts[1]

def parse_stats(path):
    # mris_anatomical_stats columns: StructName NumVert SurfArea GrayVol ThickAvg ThickStd ...
    data = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            p = line.split()
            if len(p) < 6:
                continue
            data[p[0]] = {
                'nvox': int(p[1]), 'vol_mm3': float(p[3]),
                'mean': float(p[4]), 'sd': float(p[5])
            }
    return data

def fmt(x):
    return "NaN" if x != x else f"{x:.6f}"

def mean_v(vals):
    v = [x for x in vals if x == x]
    return sum(v)/len(v) if v else float('nan')

def sd_v(vals):
    v = [x for x in vals if x == x]
    if len(v) < 2: return float('nan')
    mu = sum(v)/len(v)
    return math.sqrt(sum((x-mu)**2 for x in v)/(len(v)-1))

def z(x, mu, sd):
    return float('nan') if (x!=x or mu!=mu or sd!=sd or sd<=0) else (x-mu)/sd

lh = parse_stats(lh_stats)
rh = parse_stats(rh_stats)

lh_means, rh_means, rows = [], [], []
for idx in range(1, 49):
    name = idx2ctab.get(idx, f"ROI_{idx}")      # matches stats file StructName
    display_name = idx2display.get(idx, name)  # human-readable for output CSV
    lh_d = lh.get(name, {})
    rh_d = rh.get(name, {})
    lm = float(lh_d.get('mean', 'nan') if lh_d else 'nan')
    rm = float(rh_d.get('mean', 'nan') if rh_d else 'nan')
    lh_means.append(lm)
    rh_means.append(rm)
    rows.append({'idx': idx, 'name': name, 'display': display_name,
        'lm': lm, 'ls': float(lh_d.get('sd', 'nan') if lh_d else 'nan'),
        'ln': lh_d.get('nvox', 0), 'lv': float(lh_d.get('vol_mm3', 'nan') if lh_d else 'nan'),
        'rm': rm, 'rs': float(rh_d.get('sd', 'nan') if rh_d else 'nan'),
        'rn': rh_d.get('nvox', 0), 'rv': float(rh_d.get('vol_mm3', 'nan') if rh_d else 'nan'),
    })

all_m = lh_means + rh_means
mu_a, sd_a = mean_v(all_m),    sd_v(all_m)
mu_l, sd_l = mean_v(lh_means), sd_v(lh_means)
mu_r, sd_r = mean_v(rh_means), sd_v(rh_means)

fields = ['Index','ROIName',
    'LH_mean','LH_sd','LH_nvox','LH_vol_mm3',
    'RH_mean','RH_sd','RH_nvox','RH_vol_mm3',
    'Mean_CT','LI',
    'z_LH_global','z_RH_global','z_LH_hemi','z_RH_hemi']

out_rows = []
for r in rows:
    lm, rm = r['lm'], r['rm']
    if lm==lm and rm==rm:
        mct = (lm+rm)/2; li = (lm-rm)/(lm+rm) if (lm+rm)!=0 else float('nan')
    elif lm==lm: mct, li = lm, float('nan')
    elif rm==rm: mct, li = rm, float('nan')
    else:        mct, li = float('nan'), float('nan')
    out_rows.append({
        'Index': r['idx'], 'ROIName': r['display'],
        'LH_mean': fmt(lm), 'LH_sd': fmt(r['ls']), 'LH_nvox': r['ln'], 'LH_vol_mm3': fmt(r['lv']),
        'RH_mean': fmt(rm), 'RH_sd': fmt(r['rs']), 'RH_nvox': r['rn'], 'RH_vol_mm3': fmt(r['rv']),
        'Mean_CT': fmt(mct), 'LI': fmt(li),
        'z_LH_global': fmt(z(lm,mu_a,sd_a)), 'z_RH_global': fmt(z(rm,mu_a,sd_a)),
        'z_LH_hemi':   fmt(z(lm,mu_l,sd_l)), 'z_RH_hemi':   fmt(z(rm,mu_r,sd_r)),
    })

with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(out_rows)

n_lh = sum(1 for r in out_rows if r['LH_mean'] != 'NaN')
n_rh = sum(1 for r in out_rows if r['RH_mean'] != 'NaN')
print(f"[OK] Wrote {out_path}")
print(f"[INFO] ROIs with data — LH: {n_lh}/48  RH: {n_rh}/48")
PY

# --- Step 5: Subject volumes from aseg.stats ---------------------------------
python3 - "$STATS_DIR/aseg.stats" "$VOLCSV" <<'PY'
import sys, re, csv, os
aseg_path, out_path = sys.argv[1], sys.argv[2]
pat = re.compile(r'^\s*#\s*Measure\s+([^,]+)\s*,\s*([^,]+)\s*,[^,]*,\s*([0-9eE.+-]+)\s*,')
vals = {}
with open(aseg_path, 'r', errors='ignore') as f:
    for line in f:
        m = pat.match(line)
        if not m: continue
        a, b, v = m.group(1).strip(), m.group(2).strip(), float(m.group(3))
        vals[a] = v; vals[b] = v
def get(*keys):
    for k in keys:
        if k in vals: return vals[k]
    return None
eTIV = get('EstimatedTotalIntraCranialVol','eTIV')
GM   = get('TotalGrayVol','TotalGray')
WM   = get('CerebralWhiteMatterVol','CerebralWhiteMatter')
BSV  = get('BrainSegVol','BrainSeg')
CSF_total = (BSV-(GM+WM)) if all(x is not None for x in [BSV,GM,WM]) else None
CSF_vent  = get('VentricleChoroidVol')
def fmt(x):
    try: return f"{float(x):.6f}"
    except: return "NaN"
with open(out_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Measure','Value_mm3'])
    w.writerow(['EstimatedTotalIntraCranialVol', fmt(eTIV)])
    w.writerow(['TotalGrayVol',                  fmt(GM)])
    w.writerow(['CerebralWhiteMatterVol',         fmt(WM)])
    w.writerow(['CSF_total',                      fmt(CSF_total)])
    w.writerow(['CSF_VentricleChoroidVol',         fmt(CSF_vent)])
print("[OK] Wrote", out_path)
PY

# --- Step 6: Subject volumes vs norms ----------------------------------------
if [ ! -f "$VOLNORMS_CSV" ]; then
  echo "[SKIP] Normative volumes CSV not found — skipping volumes-vs-norms"
else
VOLCSV="$VOLCSV" VOLNORMS_CSV="$VOLNORMS_CSV" VOLPVSN_CSV="$VOLPVSN_CSV" python3 - <<'VOLS_PY'
import csv, os

vol_csv   = os.environ["VOLCSV"]
norms_csv = os.environ["VOLNORMS_CSV"]
out_csv   = os.environ["VOLPVSN_CSV"]

def tof(x):
    try:
        s = str(x).strip()
        return float("nan") if s.lower() in ("","nan","na","none") else float(s)
    except:
        return float("nan")

def zs(x, mu, sd):
    return float("nan") if (x!=x or mu!=mu or sd!=sd or sd<=0) else (x-mu)/sd

def rzs(x, med, iqr):
    if x!=x or med!=med or iqr!=iqr: return float("nan")
    sc = iqr/1.349 if iqr>0 else float("nan")
    return float("nan") if (sc!=sc or sc<=0) else (x-med)/sc

def f(x): return "NaN" if x!=x else f"{x:.6f}"

norms = {}
with open(norms_csv, newline="") as fh:
    for row in csv.DictReader(fh):
        norms[row["Measure"]] = {
            "mean":   tof(row.get("mean")),
            "sd":     tof(row.get("sd")),
            "median": tof(row.get("median")),
            "iqr":    tof(row.get("iqr")),
        }

subj = {}
with open(vol_csv, newline="") as fh:
    for row in csv.DictReader(fh):
        subj[row["Measure"]] = tof(row.get("Value_mm3"))

fields = ["Measure","Value_mm3","z","rz"]
out_rows = []
for measure in ["EstimatedTotalIntraCranialVol","TotalGrayVol",
                "CerebralWhiteMatterVol","CSF_total","CSF_VentricleChoroidVol"]:
    val = subj.get(measure, float("nan"))
    n   = norms.get(measure, {})
    out_rows.append({
        "Measure":   measure,
        "Value_mm3": f(val),
        "z":         f(zs(val, n.get("mean",float("nan")), n.get("sd",float("nan")))),
        "rz":        f(rzs(val, n.get("median",float("nan")), n.get("iqr",float("nan")))),
    })

with open(out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader(); w.writerows(out_rows)
print(f"[OK] Wrote {out_csv}")
VOLS_PY
fi

# --- Step 7: Patient vs norms ------------------------------------------------
# NOTE: NORMS_CSV must be regenerated with this surface-based pipeline before
# these comparisons are valid.
if [ ! -f "$NORMS_CSV" ]; then
  echo "[SKIP] Normative CSV not found — skipping patient-vs-norms"
else
FINAL_CSV="$FINAL_CSV" NORMS_CSV="$NORMS_CSV" PVSN_CSV="$PVSN_CSV" python3 - <<'PY'
import csv, os

patient_csv = os.environ["FINAL_CSV"]
norms_csv   = os.environ["NORMS_CSV"]
out_csv     = os.environ["PVSN_CSV"]

def tof(x):
    try:
        s = str(x).strip()
        return float('nan') if s.lower() in ('','nan','na','none') else float(s)
    except: return float('nan')

def zs(x, mu, sd):
    return float('nan') if (x!=x or mu!=mu or sd!=sd or sd<=0) else (x-mu)/sd

def rzs(x, med, iqr):
    if x!=x or med!=med or iqr!=iqr: return float('nan')
    sc = iqr/1.349 if iqr>0 else float('nan')
    return float('nan') if (sc!=sc or sc<=0) else (x-med)/sc

norms = {}
with open(norms_csv, newline='') as f:
    for row in csv.DictReader(f):
        idx = int(row["Index"])
        norms[idx] = {k: tof(row.get(k)) for k in [
            "LH_mean","LH_sd","LH_median","LH_iqr",
            "RH_mean","RH_sd","RH_median","RH_iqr",
            "Both_mean","Both_sd","Both_median","Both_iqr",
            "LI_mean","LI_sd","LI_median","LI_iqr"]}
        norms[idx]["ROIName"] = row.get("ROIName", f"ROI_{idx}")

subj = {}
with open(patient_csv, newline='') as f:
    for row in csv.DictReader(f):
        idx = int(row["Index"])
        subj[idx] = {"ROIName": row.get("ROIName",""), "LH": tof(row.get("LH_mean")), "RH": tof(row.get("RH_mean"))}

fields = ["Index","ROIName",
    "LH_patient","LH_z","LH_rz",
    "RH_patient","RH_z","RH_rz",
    "Both_patient","Both_z","Both_rz",
    "LI_patient","LI_z","LI_rz"]

def f(x): return "NaN" if x!=x else f"{x:.6f}"

out_rows = []
for idx in sorted(norms):
    n = norms[idx]
    s = subj.get(idx, {"ROIName": n["ROIName"], "LH": float('nan'), "RH": float('nan')})
    LH, RH = s["LH"], s["RH"]
    Both = (LH+RH)/2 if LH==LH and RH==RH else float('nan')
    LI   = (LH-RH)/(LH+RH) if LH==LH and RH==RH and (LH+RH)!=0 else float('nan')
    out_rows.append({
        "Index": idx, "ROIName": n["ROIName"] or s["ROIName"],
        "LH_patient": f(LH),   "LH_z": f(zs(LH,n["LH_mean"],n["LH_sd"])),     "LH_rz": f(rzs(LH,n["LH_median"],n["LH_iqr"])),
        "RH_patient": f(RH),   "RH_z": f(zs(RH,n["RH_mean"],n["RH_sd"])),     "RH_rz": f(rzs(RH,n["RH_median"],n["RH_iqr"])),
        "Both_patient": f(Both),"Both_z": f(zs(Both,n["Both_mean"],n["Both_sd"])),"Both_rz": f(rzs(Both,n["Both_median"],n["Both_iqr"])),
        "LI_patient": f(LI),   "LI_z": f(zs(LI,n["LI_mean"],n["LI_sd"])),     "LI_rz": f(rzs(LI,n["LI_median"],n["LI_iqr"])),
    })

with open(out_csv, "w", newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader(); w.writerows(out_rows)
n_lh = sum(1 for r in out_rows if r["LH_patient"]!="NaN")
n_rh = sum(1 for r in out_rows if r["RH_patient"]!="NaN")
print(f"[OK] Wrote {out_csv}")
print(f"[INFO] ROIs with data — LH: {n_lh}  RH: {n_rh}  (of {len(out_rows)})")
PY
fi

# --- Summary -----------------------------------------------------------------
echo
echo "[DONE] ${OUT_PREFIX}"
echo "  Thickness CSV:         $FINAL_CSV"
echo "  Subject volumes CSV:   $VOLCSV"
echo "  Volumes vs norms:      $VOLPVSN_CSV"
echo "  Thickness vs norms:    $PVSN_CSV"
echo "  Output folder:       $OUTDIR"
echo
echo "QC - view subject annotation in freeview:"
echo "  freeview -f $(cd "$SUBDIR" && pwd)/surf/lh.inflated:annot=${WORKDIR}/lh.ho_${THR}.annot:annot_outline=1"
