import os
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import zoom

# resize lesion masks to binary images matching atlas dimensions
def resize_binary(image, target_shape):
    if image.shape != target_shape:
        factors = [t / float(s) for s, t in zip(image.shape, target_shape)]
        resized = zoom(image, factors, order=0)
        print(f"Resized image from {image.shape} to {target_shape}")
    else:
        resized = image
    resized[resized >= 0.5] = 1
    resized[resized < 0.5] = 0
    return resized

# File paths
lesion_folder = 'path to lesion folder'  # Folder containing lesion masks (.nii files)
AAL_path = 'path to ROI files'  # Folder containing ROI images (.nii files)
output_directory_path = 'output directory'  # Output directory path

# Extract lesion files
lesion_files = [f for f in os.listdir(lesion_folder) if f.endswith('.nii') and not f.startswith('.')]
subjects = [os.path.splitext(f)[0] for f in lesion_files]

# Load ROIs
AAL_ROIs_files = [f for f in os.listdir(AAL_path) if f.endswith('.nii') and not f.startswith('.')]
first_roi = nib.load(os.path.join(AAL_path, AAL_ROIs_files[0])).get_fdata()
AAL_ROIs = np.zeros((*first_roi.shape, len(AAL_ROIs_files)), dtype=np.float32)

for i, roi_file in enumerate(AAL_ROIs_files):
    AAL_ROIs[:, :, :, i] = nib.load(os.path.join(AAL_path, roi_file)).get_fdata()

# Calculate total volume per ROI and determine dimensions
AAL_ROIs_vol = [np.sum(AAL_ROIs[:, :, :, i]) for i in range(len(AAL_ROIs_files))]
target_dimension = AAL_ROIs[:, :, :, 0].shape

# Calculate proportion of spared tissue for each subject per ROI
output_AAL = np.full((len(subjects), len(AAL_ROIs_files)), np.nan)
for i, filename in enumerate(lesion_files):
    lesionmap = nib.load(os.path.join(lesion_folder, filename)).get_fdata().astype(np.float32)
    lesionmap[np.isnan(lesionmap)] = 0
    lesionmap_resized = resize_binary(lesionmap, target_dimension)

    for h in range(len(AAL_ROIs_files)):
        overlap = lesionmap_resized * AAL_ROIs[:, :, :, h]
        vol_overlap = np.sum(overlap)
        proportion_spared = (AAL_ROIs_vol[h] - vol_overlap) / AAL_ROIs_vol[h]
        output_AAL[i, h] = proportion_spared

    print(f"Completed {subjects[i]}")

# Save results
output_AAL[np.isnan(output_AAL)] = 0
roi_names = [os.path.splitext(name)[0] for name in AAL_ROIs_files]
output_df = pd.DataFrame(output_AAL, columns=roi_names, index=subjects)
output_df.to_csv(os.path.join(output_directory_path, 'gray_matter_spared_test.csv'), index=True)
