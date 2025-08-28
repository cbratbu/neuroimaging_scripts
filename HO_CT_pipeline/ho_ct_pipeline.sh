#!/usr/bin/env bash
# HO(prob)->T1 + cortical thickness pipeline (LH/RH weighted stats + Mean_CT, LI, z-scores)
# command to run: bash /projectnb/skiran/ekropp/R01_computational_modeling/HO_CT_pipeline/ho_ct_pipeline.sh <recon-all output directory>


set -e

# --- Args --------------------------------------------------------------------
if [ $# -lt 1 ]; then
  echo "Usage: $0 /path/to/<subject_dir>   # dir with mri/, surf/, stats/"
  exit 1
fi
SUBDIR="$1"
[ -d "$SUBDIR/mri" ] || { echo "ERROR: not a FreeSurfer subject dir: $SUBDIR"; exit 1; }
echo "[INFO] Subject dir: $SUBDIR"

# --- Environment --------------------------------------------------------------
if type module >/dev/null 2>&1; then
  module load freesurfer/7.4.1 2>/dev/null || module load freesurfer 2>/dev/null || true
  module load fsl 2>/dev/null || true
fi
echo "[INFO] FREESURFER_HOME=${FREESURFER_HOME:-<unset>}"
echo "[INFO] FSLDIR=${FSLDIR:-<unset>}"
command -v flirt >/dev/null || { echo "ERROR: FLIRT not found on PATH"; exit 1; }
command -v fnirt >/dev/null || { echo "ERROR: FNIRT not found on PATH"; exit 1; }
command -v mri_surf2vol >/dev/null || { echo "ERROR: FreeSurfer binaries not on PATH"; exit 1; }

# --- Globals ------------------------------------------------------------------
GLOBAL_DIR="/projectnb/skiran/ekropp/R01_computational_modeling/HO_CT_pipeline"
LUT="${GLOBAL_DIR}/HO48_1based.lut"
LABELS_TSV="${GLOBAL_DIR}/HO48_labels_1based.tsv"
[ -f "$LUT" ] || { echo "ERROR: Missing LUT: $LUT"; exit 1; }
[ -f "$LABELS_TSV" ] || { echo "ERROR: Missing labels TSV: $LABELS_TSV"; exit 1; }
echo "[INFO] Using LUT:  $LUT"
echo "[INFO] Using TSV:  $LABELS_TSV"

# --- Subject paths ------------------------------------------------------------
export SUBJECTS_DIR="$(cd "$SUBDIR/.." && pwd)"
SUB="$(basename "$SUBDIR")"
MRI="$SUBDIR/mri"
SURF="$SUBDIR/surf"
STATS="$SUBDIR/stats"
T1MGZ="$MRI/T1.mgz"
DERIVED="$SUBDIR/derived/ho_ct"
mkdir -p "$DERIVED"

# --- Step 1: T1.nii.gz + affine + FNIRT (MNI->subject) -----------------------
T1NII="$MRI/T1.nii.gz"
[ -f "$T1NII" ] || mri_convert "$T1MGZ" "$T1NII" >/dev/null

AFF="$DERIVED/mni2sub_aff.mat"
if [ ! -f "$AFF" ]; then
  echo "[RUN] FLIRT affine (MNI->subject) -> $AFF"
  flirt -in "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" -ref "$T1NII" \
        -omat "$AFF" -dof 12 -cost mutualinfo \
        -searchrx -90 90 -searchry -90 90 -searchrz -90 90
fi

T1_BET_BASE="$DERIVED/T1_brain"
T1_BET="$T1_BET_BASE.nii.gz"
MASK="${T1_BET_BASE}_mask.nii.gz"
if [ ! -f "$MASK" ]; then
  echo "[RUN] BET for FNIRT"
  bet "$T1NII" "$T1_BET_BASE" -R -f 0.3 -g 0 -m >/dev/null
fi

NL_WARP="$DERIVED/mni2sub_warpcoef.nii.gz"
if [ ! -f "$NL_WARP" ]; then
  echo "[RUN] FNIRT nonlinear (MNI->subject) -> $NL_WARP"
  fnirt --in="$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" \
        --aff="$AFF" --cout="$NL_WARP" \
        --ref="$T1NII" --refmask="$MASK" \
        --config=T1_2_MNI152_2mm
fi

# --- Step 2: Apply HO probabilistic atlas + make HO_SUB (max-prob labels) ----
HO_PROB_MNI="$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-prob-1mm.nii.gz"
HO_PROB_SUB="$DERIVED/ho_cort_prob_sub.nii.gz"
applywarp --in="$HO_PROB_MNI" --ref="$T1NII" --warp="$NL_WARP" \
          --out="$HO_PROB_SUB" --interp=spline

HO_SUB="$DERIVED/ho_cort_thr0_sub.nii.gz"
TMP_MAXIDX="$DERIVED/ho_prob_maxidx_0based.nii.gz"
fslmaths "$HO_PROB_SUB" -Tmaxn "$TMP_MAXIDX" >/dev/null
fslmaths "$TMP_MAXIDX" -add 1 "$HO_SUB" >/dev/null
rm -f "$TMP_MAXIDX"

# --- Step 3: Hemisphere masks from ribbon (3=LH, 42=RH) ----------------------
LH_RIB_MGZ="$DERIVED/lh.ribbon.mask.mgz"
RH_RIB_MGZ="$DERIVED/rh.ribbon.mask.mgz"
LH_RIB_NII="$DERIVED/lh.ribbon.mask.nii.gz"
RH_RIB_NII="$DERIVED/rh.ribbon.mask.nii.gz"
mri_binarize --i "$MRI/ribbon.mgz" --match 3  --o "$LH_RIB_MGZ"
mri_binarize --i "$MRI/ribbon.mgz" --match 42 --o "$RH_RIB_MGZ"
mri_convert "$LH_RIB_MGZ" "$LH_RIB_NII" >/dev/null
mri_convert "$RH_RIB_MGZ" "$RH_RIB_NII" >/dev/null

# --- Step 4: Rasterize surface thickness to T1 (identity LTA) ----------------
LTA="$DERIVED/surf2vol.identity.lta"
lta_convert --inlta identity.nofile --outlta "$LTA" --src "$T1MGZ" --trg "$T1MGZ"
LH_THICK_MGZ="$DERIVED/lh.thickness.ribbon.mgz"
RH_THICK_MGZ="$DERIVED/rh.thickness.ribbon.mgz"
LH_THICK_NII="$DERIVED/lh.thickness.ribbon.nii.gz"
RH_THICK_NII="$DERIVED/rh.thickness.ribbon.nii.gz"

if [ ! -f "$LH_THICK_MGZ" ]; then
  echo "[RUN] mri_surf2vol (LH thickness)"
  mri_surf2vol --subject "$SUB" --hemi lh \
               --surfval "$SURF/lh.thickness" \
               --reg "$LTA" --template "$T1MGZ" \
               --fillribbon --o "$LH_THICK_MGZ"
fi
if [ ! -f "$RH_THICK_MGZ" ]; then
  echo "[RUN] mri_surf2vol (RH thickness)"
  mri_surf2vol --subject "$SUB" --hemi rh \
               --surfval "$SURF/rh.thickness" \
               --reg "$LTA" --template "$T1MGZ" \
               --fillribbon --o "$RH_THICK_MGZ"
fi
mri_convert "$LH_THICK_MGZ" "$LH_THICK_NII" >/dev/null
mri_convert "$RH_THICK_MGZ" "$RH_THICK_NII" >/dev/null

# --- Step 5: Probability-weighted ROI stats per hemi -------------------------
RAW_CSV="$DERIVED/HO48_thickness_means_raw.csv"
echo "Index,ROIName,LH_mean,LH_sd_w,LH_min,LH_max,LH_nvox,LH_vol_mm3,RH_mean,RH_sd_w,RH_min,RH_max,RH_nvox,RH_vol_mm3" > "$RAW_CSV"

dx=$(fslval "$T1NII" pixdim1); dy=$(fslval "$T1NII" pixdim2); dz=$(fslval "$T1NII" pixdim3)
voxvol=$(python3 - <<PY
dx=float("$dx"); dy=float("$dy"); dz=float("$dz")
print(f"{dx*dy*dz:.6f}")
PY
)

LH_THICK2_NII="$DERIVED/lh.thickness2.ribbon.nii.gz"
RH_THICK2_NII="$DERIVED/rh.thickness2.ribbon.nii.gz"
fslmaths "$LH_THICK_NII" -sqr "$LH_THICK2_NII" >/dev/null
fslmaths "$RH_THICK_NII" -sqr "$RH_THICK2_NII" >/dev/null

for IDX in $(seq 1 48); do
  NAME=$(awk -v i="$IDX" -F'\t' '$1==i {print $2}' "$LABELS_TSV")

  ROI_PROB="$DERIVED/roiProb_${IDX}.nii.gz"
  fslroi "$HO_PROB_SUB" "$ROI_PROB" $((IDX-1)) 1

  ROI_L="$DERIVED/roiProb_${IDX}_LH.nii.gz"
  ROI_R="$DERIVED/roiProb_${IDX}_RH.nii.gz"
  fslmaths "$ROI_PROB" -mul "$LH_RIB_NII" "$ROI_L" >/dev/null
  fslmaths "$ROI_PROB" -mul "$RH_RIB_NII" "$ROI_R" >/dev/null

  ROI_L_BIN="$DERIVED/roiProb_${IDX}_LH_bin.nii.gz"
  ROI_R_BIN="$DERIVED/roiProb_${IDX}_RH_bin.nii.gz"
  fslmaths "$ROI_L" -thr 0 -bin "$ROI_L_BIN" >/dev/null
  fslmaths "$ROI_R" -thr 0 -bin "$ROI_R_BIN" >/dev/null

  NUM_L="$DERIVED/tmp_num_L_${IDX}.nii.gz"; DEN_L="$DERIVED/tmp_den_L_${IDX}.nii.gz"; NUM2_L="$DERIVED/tmp_num2_L_${IDX}.nii.gz"
  fslmaths "$LH_THICK_NII"  -mul "$ROI_L" "$NUM_L"  >/dev/null
  fslmaths "$LH_THICK2_NII" -mul "$ROI_L" "$NUM2_L" >/dev/null
  fslmaths "$ROI_L" -mul 1 "$DEN_L" >/dev/null
  LH_NUM=$(fslstats "$NUM_L"  -M -V | awk '{print $1*$3}')
  LH_NUM2=$(fslstats "$NUM2_L" -M -V | awk '{print $1*$3}')
  LH_DEN=$(fslstats "$DEN_L"  -M -V | awk '{print $1*$3}')
  LH_MEAN=$(python3 - <<PY
num=float("$LH_NUM"); den=float("$LH_DEN")
print("NaN" if den<=0 else f"{num/den:.6f}")
PY
)
  LH_SD_W=$(python3 - <<PY
num=float("$LH_NUM"); num2=float("$LH_NUM2"); den=float("$LH_DEN")
import math
if den<=0: print("NaN")
else:
    mean = num/den
    var = num2/den - mean*mean
    if var < 0: var = 0.0
    print(f"{math.sqrt(var):.6f}")
PY
)
  LH_MINMAX=$(fslstats "$LH_THICK_NII" -k "$ROI_L_BIN" -R 2>/dev/null || echo "nan nan")
  LH_MIN=$(echo "$LH_MINMAX" | awk '{print $1}')
  LH_MAX=$(echo "$LH_MINMAX" | awk '{print $2}')
  LH_NVOX=$(fslstats "$ROI_L_BIN" -V | awk '{print $1}')
  LH_VMM3=$(python3 - <<PY
den=float("$LH_DEN"); vv=float("$voxvol")
print("NaN" if den<=0 else f"{den*vv:.6f}")
PY
)

  NUM_R="$DERIVED/tmp_num_R_${IDX}.nii.gz"; DEN_R="$DERIVED/tmp_den_R_${IDX}.nii.gz"; NUM2_R="$DERIVED/tmp_num2_R_${IDX}.nii.gz"
  fslmaths "$RH_THICK_NII"  -mul "$ROI_R" "$NUM_R"  >/dev/null
  fslmaths "$RH_THICK2_NII" -mul "$ROI_R" "$NUM2_R" >/dev/null
  fslmaths "$ROI_R" -mul 1 "$DEN_R" >/dev/null
  RH_NUM=$(fslstats "$NUM_R"  -M -V | awk '{print $1*$3}')
  RH_NUM2=$(fslstats "$NUM2_R" -M -V | awk '{print $1*$3}')
  RH_DEN=$(fslstats "$DEN_R"  -M -V | awk '{print $1*$3}')
  RH_MEAN=$(python3 - <<PY
num=float("$RH_NUM"); den=float("$RH_DEN")
print("NaN" if den<=0 else f"{num/den:.6f}")
PY
)
  RH_SD_W=$(python3 - <<PY
num=float("$RH_NUM"); num2=float("$RH_NUM2"); den=float("$RH_DEN")
import math
if den<=0: print("NaN")
else:
    mean = num/den
    var = num2/den - mean*mean
    if var < 0: var = 0.0
    print(f"{math.sqrt(var):.6f}")
PY
)
  RH_MINMAX=$(fslstats "$RH_THICK_NII" -k "$ROI_R_BIN" -R 2>/dev/null || echo "nan nan")
  RH_MIN=$(echo "$RH_MINMAX" | awk '{print $1}')
  RH_MAX=$(echo "$RH_MINMAX" | awk '{print $2}')
  RH_NVOX=$(fslstats "$ROI_R_BIN" -V | awk '{print $1}')
  RH_VMM3=$(python3 - <<PY
den=float("$RH_DEN"); vv=float("$voxvol")
print("NaN" if den<=0 else f"{den*vv:.6f}")
PY
)

  echo "${IDX},\"${NAME}\",${LH_MEAN},${LH_SD_W},${LH_MIN},${LH_MAX},${LH_NVOX},${LH_VMM3},${RH_MEAN},${RH_SD_W},${RH_MIN},${RH_MAX},${RH_NVOX},${RH_VMM3}" >> "$RAW_CSV"

  rm -f "$ROI_PROB" "$ROI_L" "$ROI_R" "$ROI_L_BIN" "$ROI_R_BIN" "$NUM_L" "$DEN_L" "$NUM2_L" "$NUM_R" "$DEN_R" "$NUM2_R"
done

# --- Step 6: metrics (Mean_CT, LI, z-scores) ---------------------------------
FINAL_CSV="$DERIVED/HO48_thickness_metrics.csv"
DERIVED="$DERIVED" python3 - <<'PY'
import csv, os, math
der=os.environ['DERIVED']
infile=os.path.join(der,'HO48_thickness_means_raw.csv')
outfile=os.path.join(der,'HO48_thickness_metrics.csv')
def to_f(x):
    try:
        s=str(x).strip().lower()
        if s in ('','nan','na','none'): return float('nan')
        return float(s)
    except: return float('nan')
def mean(a):
    d=[v for v in a if v==v]; return sum(d)/len(d) if d else float('nan')
def sd(a):
    d=[v for v in a if v==v]; n=len(d)
    if n<2: return float('nan')
    mu=sum(d)/n
    return (sum((v-mu)**2 for v in d)/(n-1))**0.5
rows=list(csv.DictReader(open(infile)))
LH=[to_f(r['LH_mean']) for r in rows]
RH=[to_f(r['RH_mean']) for r in rows]
mu_all, sd_all = mean(LH+RH), sd(LH+RH)
mu_L,   sd_L   = mean(LH),    sd(LH)
mu_R,   sd_R   = mean(RH),    sd(RH)
out=[]
for r in rows:
    lh, rh = to_f(r['LH_mean']), to_f(r['RH_mean'])
    if lh==lh and rh==rh: m=(lh+rh)/2.0
    elif lh==lh:          m=lh
    elif rh==rh:          m=rh
    else:                 m=float('nan')
    li = (lh-rh)/(lh+rh) if (lh==lh and rh==rh and (lh+rh)!=0) else float('nan')
    out.append({
        'Index'        : r['Index'],
        'ROIName'      : r['ROIName'],
        'LH_mean'      : 'NaN' if lh!=lh else f'{lh:.6f}',
        'RH_mean'      : 'NaN' if rh!=rh else f'{rh:.6f}',
        'Mean_CT'      : 'NaN' if m!=m  else f'{m:.6f}',
        'LI'           : 'NaN' if li!=li else f'{li:.6f}',
        'z_LH_global'  : 'NaN' if (lh!=lh or sd_all!=sd_all or sd_all==0) else f'{(lh-mu_all)/sd_all:.6f}',
        'z_RH_global'  : 'NaN' if (rh!=rh or sd_all!=sd_all or sd_all==0) else f'{(rh-mu_all)/sd_all:.6f}',
        'z_LH_hemi'    : 'NaN' if (lh!=lh or sd_L!=sd_L   or sd_L==0)   else f'{(lh-mu_L)/sd_L:.6f}',
        'z_RH_hemi'    : 'NaN' if (rh!=rh or sd_R!=sd_R   or sd_R==0)   else f'{(rh-mu_R)/sd_R:.6f}',
    })
with open(outfile,'w',newline='') as f:
    w=csv.DictWriter(f, fieldnames=['Index','ROIName','LH_mean','RH_mean','Mean_CT','LI','z_LH_global','z_RH_global','z_LH_hemi','z_RH_hemi'])
    w.writeheader(); w.writerows(out)
print("[OK] Wrote", outfile)
PY

# --- Step 7: subject volumes (eTIV, GM, WM, CSF_total, CSF_VentricleChoroid) --
VOLCSV="$DERIVED/subject_volumes.csv"
python3 - "$STATS/aseg.stats" "$VOLCSV" <<'PY'
import sys, re, csv, os
aseg_path, out_path = sys.argv[1], sys.argv[2]
pat = re.compile(r'^\s*#\s*Measure\s+([^,]+)\s*,\s*([^,]+)\s*,[^,]*,\s*([0-9eE.+-]+)\s*,')
vals={}
with open(aseg_path,'r',errors='ignore') as f:
    for line in f:
        m=pat.match(line)
        if not m: continue
        a,b,v=m.group(1).strip(),m.group(2).strip(),float(m.group(3))
        vals[a]=v; vals[b]=v
def get(*keys):
    for k in keys:
        if k in vals: return vals[k]
    return None
eTIV = get('EstimatedTotalIntraCranialVol','eTIV')
GM   = get('TotalGrayVol','TotalGray')
WM   = get('CerebralWhiteMatterVol','CerebralWhiteMatter')
BSV  = get('BrainSegVol','BrainSeg')
CSF_total = (BSV - (GM + WM)) if (BSV is not None and GM is not None and WM is not None) else None
CSF_vent  = get('VentricleChoroidVol')
def fmt(x):
    try: return f"{float(x):.6f}"
    except: return "NaN"
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['Measure','Value_mm3'])
    w.writerow(['EstimatedTotalIntraCranialVol', fmt(eTIV)])
    w.writerow(['TotalGrayVol', fmt(GM)])
    w.writerow(['CerebralWhiteMatterVol', fmt(WM)])
    w.writerow(['CSF_total', fmt(CSF_total)])
    w.writerow(['CSF_VentricleChoroidVol', fmt(CSF_vent)])
print("[OK] Wrote", out_path)
PY

# --- Summary & Freeview (cortex-masked labels) --------------------------------
# Build a cortex-only mask and a masked label image for visualization
CORTEX_MASK="$DERIVED/cortex_mask.nii.gz"
HO_LABELS_CTX="$DERIVED/ho_labels_cortex_only.nii.gz"
fslmaths "$DERIVED/lh.ribbon.mask.nii.gz" -add "$DERIVED/rh.ribbon.mask.nii.gz" -bin "$CORTEX_MASK" >/dev/null
fslmaths "$HO_SUB" -mas "$CORTEX_MASK" "$HO_LABELS_CTX" >/dev/null

echo
echo " Finished. Outputs in: $DERIVED"
echo " - Atlas label (for Freeview):   $HO_SUB"
echo " - Prob atlas (for weighting):   $HO_PROB_SUB"
echo " - LH/RH thickness NIfTI:        $LH_THICK_NII , $RH_THICK_NII"
echo " - ROI stats CSV:                $RAW_CSV"
echo " - Metrics CSV:                  $FINAL_CSV"
echo " - Subject volumes CSV:          $VOLCSV"
echo
echo "Freeview QC (cortex-only labels):"
echo "freeview -v $MRI/T1.mgz \\"
echo "  $HO_LABELS_CTX:colormap=lut:lut=$LUT:opacity=0.45 \\"
echo "  $MRI/brainmask.mgz:opacity=0.15"