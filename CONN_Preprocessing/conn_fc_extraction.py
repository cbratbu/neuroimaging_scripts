#!/usr/bin/env python3

################
# Usage: python3 extract_conn_fc.py /path/to/conn_fc_export.nii[.gz] (with .json file of matching name from conn export saved to same location)
# Output: /fc_outputs_<nifti_basename>/ containing one CSV per seed ROI
# Each CSV: rows = subjects, columns = target ROIs (all ROIs except the seed for that file)
################


import os
import sys
import json
import re
import pandas as pd
import nibabel as nib

def usage_exit():
    prog = os.path.basename(sys.argv[0])
    print(f"Usage: {prog} /full/path/to/conn_fc.nii", file=sys.stderr)
    print("Assumes a JSON sidecar with the same basename: conn_fc.json", file=sys.stderr)
    sys.exit(1)

def load_conn_cube(nii_path: str, json_path: str):
    img = nib.load(nii_path)
    arr = img.get_fdata()

    # Expect (nROI, nROI, 1, nSubj) or (nROI, nROI, nSubj)
    if arr.ndim == 4 and arr.shape[2] == 1:
        arr = arr[:, :, 0, :]
    elif arr.ndim == 3:
        pass
    else:
        raise ValueError(
            f"Unexpected NIfTI shape {arr.shape}; "
            "expected (nROI, nROI, 1, nSubj) or (nROI, nROI, nSubj)"
        )

    with open(json_path, "r") as f:
        meta = json.load(f)

    roi_names = meta["names"]
    subjects = [s.replace(" rest", "").strip() for s in meta["samples"]]
    return arr, roi_names, subjects

def clean_label(name: str) -> str:
    name = re.sub(r"^atlas\.", "", name)
    name = re.sub(r"\s*\(.*?\)\s*$", "", name)
    return name

def safe_filename(s: str) -> str:
    s = s.strip()
    s = re.sub(r"[\/\\:\*\?\"<>\|]+", "_", s)
    s = re.sub(r"\s+", "_", s)
    return s

def extract_seed_df(cube, roi_names, subjects, seed_idx: int) -> pd.DataFrame:
    mat = cube[seed_idx, :, :].T  # (nSubj, nROI)

    keep_cols = [i for i in range(len(roi_names)) if i != seed_idx]
    mat = mat[:, keep_cols]
    col_names = [clean_label(roi_names[i]) for i in keep_cols]

    df = pd.DataFrame(mat, columns=col_names)
    df.insert(0, "Subject", subjects)
    return df

def main():
    if len(sys.argv) != 2:
        usage_exit()

    nii_path = sys.argv[1]
    if not os.path.isfile(nii_path):
        raise FileNotFoundError(f"NIfTI not found: {nii_path}")

    base, ext = os.path.splitext(nii_path)
    if ext.lower() == ".gz" and base.lower().endswith(".nii"):
        base = os.path.splitext(base)[0]  # strip .nii from .nii.gz
    json_path = base + ".json"

    if not os.path.isfile(json_path):
        raise FileNotFoundError(
            f"JSON sidecar not found: {json_path}\n"
            "Expected same basename as NIfTI with .json extension."
        )

    nii_dir = os.path.dirname(os.path.abspath(nii_path))
    nii_stem = os.path.basename(base)
    out_dir = os.path.join(nii_dir, f"fc_outputs_{safe_filename(nii_stem)}")
    os.makedirs(out_dir, exist_ok=True)

    cube, roi_names, subjects = load_conn_cube(nii_path, json_path)

    if cube.shape[0] != cube.shape[1]:
        raise ValueError(f"ROIxROI matrix should be square; got {cube.shape[:2]}")
    if cube.shape[0] != len(roi_names):
        raise ValueError(
            f"Mismatch: cube has {cube.shape[0]} ROIs but JSON has {len(roi_names)} names"
        )

    print(f"Loaded: {nii_path}")
    print(f"JSON:    {json_path}")
    print(f"Cube shape: {cube.shape} (nROI, nROI, nSubj)")
    print(f"ROIs: {len(roi_names)} | Subjects: {len(subjects)}")
    print(f"Writing seed CSVs to: {out_dir}")

    for seed_idx, seed_fullname in enumerate(roi_names):
        seed_label = clean_label(seed_fullname)
        df_seed = extract_seed_df(cube, roi_names, subjects, seed_idx)

        out_name = f"FC_seed_{seed_idx:03d}_{safe_filename(seed_label)}.csv"
        out_path = os.path.join(out_dir, out_name)
        df_seed.to_csv(out_path, index=False)

        if (seed_idx + 1) % 25 == 0 or seed_idx == 0 or seed_idx == len(roi_names) - 1:
            print(f"[OK] {seed_idx+1}/{len(roi_names)} saved: {out_name} | shape: {df_seed.shape}")

    print("Done.")

if __name__ == "__main__":
    main()
