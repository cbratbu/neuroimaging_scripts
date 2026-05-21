#!/usr/bin/env python3

################
# Usage:
#   python3 lesion_load.py /path/to/lesion.nii[.gz] 
#
# *Use "out=/DIR" to specify output directory, otherwise CSV is saved next to lesion file as <lesion_name>_<atlas_name>_lesionLoad.csv
#   Ex: python3 lesion_load.py /path/to/lesion.nii[.gz] out=/OUTPUT_DIRECTORY
#
# Edit ATLAS_DIR below to point to a folder containing one ROI mask per file (.nii/.nii.gz) in MNI space.
# uses first image in ATLAS_DIR as reference, so make sure all ROI images are in same space
################


import os, sys, glob
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.image import resample_to_img


ATLAS_DIR = "/projectnb/openneuro-aphasia/ML_aphasia_scripts/harvardoxford"

def strip_nii(fname: str) -> str:
    return fname.replace(".nii.gz", "").replace(".nii", "")

def finite_copy(img):
    data = img.get_fdata().astype(np.float32, copy=True)
    np.nan_to_num(data, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    return nib.Nifti1Image(data, img.affine, img.header)

def half_range_threshold_from_finite(img):
    arr = img.get_fdata()
    finite = np.isfinite(arr)
    if not np.any(finite):
        return 0.5
    vmin = float(np.min(arr[finite]))
    vmax = float(np.max(arr[finite]))
    return 0.5 if vmax == vmin else 0.5 * (vmin + vmax)

def main():

    if len(sys.argv) < 2:
        prog = os.path.basename(sys.argv[0])
        print(f"Usage: {prog} /path/to/lesion.nii[.gz] [out=/outputdir]", file=sys.stderr)
        sys.exit(1)

    lesion_path = sys.argv[1]
    if not os.path.isfile(lesion_path):
        raise FileNotFoundError(f"Lesion file not found: {lesion_path}")

    # Optional output directory
    out_dir = None
    for arg in sys.argv[2:]:
        if arg.startswith("out="):
            out_dir = arg.split("=", 1)[1]

    if not os.path.isdir(ATLAS_DIR):
        raise FileNotFoundError(f"ATLAS_DIR not found: {ATLAS_DIR}")

    roi_files = sorted(glob.glob(os.path.join(ATLAS_DIR, "*.nii*")))
    roi_files = [p for p in roi_files if os.path.isfile(p)]
    if not roi_files:
        raise FileNotFoundError(f"No ROI NIfTI files found in {ATLAS_DIR}")

    lesion_img_orig = nib.load(lesion_path)
    thr = half_range_threshold_from_finite(lesion_img_orig)
    lesion_img = finite_copy(lesion_img_orig)

    ref_img = nib.load(roi_files[0])

    if lesion_img.shape != ref_img.shape or not np.allclose(lesion_img.affine, ref_img.affine, atol=1e-5):
        lesion_img = resample_to_img(
            lesion_img, ref_img,
            interpolation="nearest",
            force_resample=True,
            copy_header=True,
            fill_value=0.0,
        )

    lesion_bin = lesion_img.get_fdata() > thr

    rows = []
    for rf in roi_files:
        roi_img = nib.load(rf)
        if roi_img.shape != ref_img.shape or not np.allclose(roi_img.affine, ref_img.affine, atol=1e-5):
            roi_img = resample_to_img(
                roi_img, ref_img,
                interpolation="nearest",
                force_resample=True,
                copy_header=True,
                fill_value=0.0,
            )

        roi_mask = roi_img.get_fdata() > 0
        roi_vox = int(np.sum(roi_mask))
        overlap = int(np.sum(lesion_bin & roi_mask)) if roi_vox else 0
        load = (overlap / roi_vox) if roi_vox else 0.0

        rows.append({"ROI": strip_nii(os.path.basename(rf)), "LesionLoad": load})

    df = pd.DataFrame(rows).sort_values("ROI")

    lesion_name = strip_nii(os.path.basename(lesion_path))
    atlas_name = os.path.basename(os.path.normpath(ATLAS_DIR))

    if out_dir is None:
        out_dir = os.path.dirname(os.path.abspath(lesion_path))

    os.makedirs(out_dir, exist_ok=True)

    out_csv = os.path.join(out_dir, f"{lesion_name}_{atlas_name}_lesionLoad.csv")

    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
