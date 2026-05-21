#!/usr/bin/env python3
# hcp_normative_ct_aggregate.py
# Aggregates per-subject HCP thickness CSVs into normative statistics CSV.
# Computes per-ROI mean, SD, median, IQR, and Shapiro-Wilk p-value across
# subjects for LH, RH, Both (mean of LH+RH), and LI ((LH-RH)/(LH+RH)).
#
# Usage:
#   python3 hcp_normative_ct_aggregate.py <output_base_dir> <norms_out_csv>

import sys
import os
import csv
import math
import glob
from collections import defaultdict

def usage():
    print("Usage: python3 hcp_normative_ct_aggregate.py <output_base_dir> <norms_out_csv>")
    sys.exit(1)

if len(sys.argv) < 3:
    usage()

OUT_BASE   = sys.argv[1]
NORMS_OUT  = sys.argv[2]
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LABELS_TSV = os.path.join(SCRIPT_DIR, "HO48_labels_1based.tsv")

# --- Load ROI name mapping ---------------------------------------------------
idx2display = {}
with open(LABELS_TSV) as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) >= 2:
            idx2display[int(p[0])] = p[1]

# --- Collect per-subject data ------------------------------------------------
# data[idx] = {'LH': [...], 'RH': [...]}
data = defaultdict(lambda: {'LH': [], 'RH': []})

csv_files = sorted(glob.glob(os.path.join(OUT_BASE, "*", "*_HO48_thickness.csv")))
n_subjects = 0

for fpath in csv_files:
    with open(fpath, newline='') as f:
        reader = csv.DictReader(f)
        sub_rows = list(reader)

    if not sub_rows:
        continue

    # Validate expected columns
    if 'LH_mean' not in sub_rows[0] or 'RH_mean' not in sub_rows[0]:
        print(f"WARNING: Unexpected format in {fpath}, skipping")
        continue

    n_subjects += 1
    for row in sub_rows:
        idx = int(row['Index'])
        try:
            lm = float(row['LH_mean'])
            if lm == lm:  # not NaN
                data[idx]['LH'].append(lm)
        except (ValueError, KeyError):
            pass
        try:
            rm = float(row['RH_mean'])
            if rm == rm:
                data[idx]['RH'].append(rm)
        except (ValueError, KeyError):
            pass

print(f"[INFO] Loaded {n_subjects} subjects")

# --- Statistical functions ---------------------------------------------------
def mean(v):
    return sum(v) / len(v) if v else float('nan')

def sd(v):
    if len(v) < 2: return float('nan')
    mu = sum(v) / len(v)
    return math.sqrt(sum((x - mu)**2 for x in v) / (len(v) - 1))

def median(v):
    if not v: return float('nan')
    s = sorted(v)
    n = len(s)
    if n % 2 == 1:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2.0

def iqr(v):
    if len(v) < 4: return float('nan')
    s = sorted(v)
    n = len(s)
    q1 = median(s[:n // 2])
    q3 = median(s[(n + 1) // 2:])
    return q3 - q1

def shapiro_wilk_p(v):
    # scipy is the standard tool for this; fall back to NaN if unavailable
    if len(v) < 3: return float('nan')
    try:
        from scipy import stats
        _, p = stats.shapiro(v)
        return float(p)
    except ImportError:
        return float('nan')

def fmt(x):
    return "NaN" if (isinstance(x, float) and x != x) else f"{x:.6f}"

# --- Build normative table ---------------------------------------------------
fields = [
    'Index', 'ROIName', 'N',
    'LH_mean', 'LH_sd', 'LH_median', 'LH_iqr', 'LH_normal',
    'RH_mean', 'RH_sd', 'RH_median', 'RH_iqr', 'RH_normal',
    'Both_mean', 'Both_sd', 'Both_median', 'Both_iqr', 'Both_normal',
    'LI_mean', 'LI_sd', 'LI_median', 'LI_iqr',
]

out_rows = []
for idx in range(1, 49):
    lh_vals = data[idx]['LH']
    rh_vals = data[idx]['RH']

    out_rows.append({
        'Index':   idx,
        'ROIName': idx2display.get(idx, f"ROI_{idx}"),
        'N':       n_subjects,
        'LH_mean':   fmt(mean(lh_vals)),   'LH_sd':   fmt(sd(lh_vals)),
        'LH_median': fmt(median(lh_vals)), 'LH_iqr':  fmt(iqr(lh_vals)),
        'LH_normal': fmt(shapiro_wilk_p(lh_vals)),
        'RH_mean':   fmt(mean(rh_vals)),   'RH_sd':   fmt(sd(rh_vals)),
        'RH_median': fmt(median(rh_vals)), 'RH_iqr':  fmt(iqr(rh_vals)),
        'RH_normal': fmt(shapiro_wilk_p(rh_vals)),
        # Both and LI require matched LH+RH pairs per subject; computed below
        'Both_mean': 'NaN', 'Both_sd': 'NaN', 'Both_median': 'NaN',
        'Both_iqr': 'NaN', 'Both_normal': 'NaN',
        'LI_mean': 'NaN', 'LI_sd': 'NaN', 'LI_median': 'NaN', 'LI_iqr': 'NaN',
    })

# Re-collect paired Both and LI values per ROI from CSVs to get matched
# LH+RH values per subject (required for correct bilateral mean and LI stats)
paired = defaultdict(lambda: {'both': [], 'li': []})
for fpath in csv_files:
    with open(fpath, newline='') as f:
        for row in csv.DictReader(f):
            idx = int(row['Index'])
            try:
                lm = float(row['LH_mean'])
                rm = float(row['RH_mean'])
                if lm == lm and rm == rm:
                    paired[idx]['both'].append((lm + rm) / 2.0)
                    denom = lm + rm
                    if denom != 0:
                        paired[idx]['li'].append((lm - rm) / denom)
            except (ValueError, KeyError):
                pass

for row in out_rows:
    idx = row['Index']
    bv = paired[idx]['both']
    lv = paired[idx]['li']
    row['Both_mean']   = fmt(mean(bv))
    row['Both_sd']     = fmt(sd(bv))
    row['Both_median'] = fmt(median(bv))
    row['Both_iqr']    = fmt(iqr(bv))
    row['Both_normal'] = fmt(shapiro_wilk_p(bv))
    row['LI_mean']     = fmt(mean(lv))
    row['LI_sd']       = fmt(sd(lv))
    row['LI_median']   = fmt(median(lv))
    row['LI_iqr']      = fmt(iqr(lv))

# --- Write output ------------------------------------------------------------
os.makedirs(os.path.dirname(os.path.abspath(NORMS_OUT)), exist_ok=True)
with open(NORMS_OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(out_rows)

print(f"[OK] Wrote {NORMS_OUT} ({len(out_rows)} ROIs, {n_subjects} subjects)")
