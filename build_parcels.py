import pandas as pd
import os

# Input paths
excel_path = "/Volumes/kiran/KiranLab4/R01 computational bilingual /1. Preparation materials.updates/Imaging materials/ROIs_per_deficit.xlsx"
roi_dir = "/Volumes/kiran/KiranLab4/R01 computational bilingual /1. Preparation materials.updates/Imaging materials/harvardoxford"
output_dir = "/Volumes/kiran/KiranLab4/R01 computational bilingual /1. Preparation materials.updates/Imaging materials/psa_ppa_impairment_ROIs"

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Read Excel
df = pd.read_excel(excel_path)

# Keep only valid ROI filenames
df = df[df['Harvard-Oxford_ROI'].astype(str).str.endswith('.nii.gz')]

# Count distinct sources per (Impairment, ROI)
counts = (
    df
    .groupby(['Impairment', 'Harvard-Oxford_ROI'])['Source']
    .nunique()
    .reset_index(name='n_sources')
)

# Filter to ROIs with ≥2 distinct sources
valid_rois = counts[counts['n_sources'] >= 2]

# Group by impairment to get sorted list of qualifying ROIs
impairment_groups = (
    valid_rois
    .groupby('Impairment')['Harvard-Oxford_ROI']
    .apply(lambda x: sorted(x))
)

# Write out shell script
with open("run_fslmaths.sh", "w") as f:
    f.write("#!/bin/bash\n\n")
    for impairment, roi_list in impairment_groups.items():
        if not roi_list:
            continue

        out_file = os.path.join(output_dir, f"{impairment}.nii.gz")
        roi_paths = [os.path.join(roi_dir, roi) for roi in roi_list]

        if len(roi_paths) == 1:
            cmd = f"fslmaths '{roi_paths[0]}' '{out_file}'"
        else:
            cmd = f"fslmaths '{roi_paths[0]}'"
            for p in roi_paths[1:]:
                cmd += f" -add '{p}'"
            cmd += f" -bin '{out_file}'"

        f.write(cmd + "\n")

print("Shell script 'run_fslmaths.sh' created. Run it with: bash run_fslmaths.sh")
