=======================================
DWI PROCESSING PIPELINE
=======================================

This folder contains scripts for preprocessing diffusion-weighted imaging (DWI) data.
The pipeline has two stages: (1) FSL preprocessing, and (2) DSI Studio processing.
Both single-subject and batch processing are supported.

=======================================
I. PIPELINE OVERVIEW
=======================================

Stage 1 — FSL preprocessing (dwi_fsl_preprocessing.sh):
  Runs TOPUP + EDDY correction on raw AP/PA DWI data.
  Outputs eddy-corrected images and a brain mask.

Stage 2 — DSI Studio pipeline (dwi_dsi_pipeline.sh):
  Runs the full DSI Studio workflow in sequence:
    Step 1: Generate SRC file (action=src)
    Step 2: GQI reconstruction (action=rec, method=4)
    Step 3: AutoTrack tractography (action=atk)
    Step 4: Aggregate per-tract statistics into a single TSV

=======================================
II. SINGLE-SUBJECT PROCESSING
=======================================

---------------------------------------
1. Required directory structure
---------------------------------------

Each subject must have a dwi/ folder containing the raw input files. Two
structures are supported: subjects with a single timepoint have dwi/ directly
under 2.preprocessing/; subjects with multiple timepoints have dwi/ under a
labeled session folder within 2.preprocessing/.

Single timepoint:
  .../sub-XXXX/2.preprocessing/dwi/

Multiple timepoints:
  .../sub-XXXX/2.preprocessing/ses-XXXX/dwi/

Required files in dwi/:

  AP direction:
    sub-XXXX_dir-AP_dwi.nii(.gz)
    sub-XXXX_dir-AP_dwi.bval
    sub-XXXX_dir-AP_dwi.bvec
    sub-XXXX_dir-AP_dwi.json

  PA direction:
    sub-XXXX_dir-PA_dwi.nii(.gz)
    sub-XXXX_dir-PA_dwi.bval
    sub-XXXX_dir-PA_dwi.bvec
    sub-XXXX_dir-PA_dwi.json

For session subjects, the session label is included in the filename:
  sub-XXXX_ses-XXXX_dir-AP_dwi.nii(.gz)  (and corresponding files)

Both .nii and .nii.gz inputs are accepted.

---------------------------------------
2. Step 1: FSL preprocessing
---------------------------------------

Script: dwi_fsl_preprocessing.sh

Usage:
  dwi_fsl_preprocessing.sh /path/to/sub-XXXX/2.preprocessing/dwi/
  dwi_fsl_preprocessing.sh /path/to/sub-XXXX/2.preprocessing/ses-XXXX/dwi/

This script performs the full FSL DWI preprocessing pipeline:
  - Reads AP/PA JSON metadata
  - Creates acqparams.txt and index.txt
  - Extracts b0 volumes from AP and PA runs
  - Runs TOPUP
  - Computes mean b0 (hifi_nodif.nii.gz) and brain mask (hifi_nodif_brain_mask.nii.gz)
  - Runs EDDY (uses eddy_cuda8.0 if GPU available, otherwise CPU eddy)

Key outputs:
  eddy_unwarped_images.nii.gz
  eddy_unwarped_images.eddy_rotated_bvecs
  hifi_nodif.nii.gz
  hifi_nodif_brain_mask.nii.gz
  topup_AP_PA_b0_* files
  eddy QC reports

These outputs are required inputs for the DSI Studio pipeline.

---------------------------------------
3. Step 2: DSI Studio pipeline
---------------------------------------

Script: dwi_dsi_pipeline.sh

Usage:
  dwi_dsi_pipeline.sh /path/to/sub-XXXX/2.preprocessing/dwi/
  dwi_dsi_pipeline.sh /path/to/sub-XXXX/2.preprocessing/ses-XXXX/dwi/

Required inputs in dwi_dir:
  eddy_unwarped_images.nii.gz              (from FSL step)
  eddy_unwarped_images.eddy_rotated_bvecs  (from FSL step)
  sub-*_dir-AP_dwi.bval
  hifi_nodif_brain_mask.nii.gz             (optional; auto-mask used if absent)

The script runs in four steps:

  Step 1 — SRC generation:
    Creates a DSI Studio source file (.sz) from the eddy-corrected NIfTI.
    Skipped if sub-XXXX_eddy_unwarped_images.sz already exists.

  Step 2 — GQI reconstruction:
    Reconstructs fiber orientation distribution functions in native diffusion space
    using Generalized Q-Sampling Imaging (GQI, method=4). GQI is a model-free method
    that resolves crossing fibers and computes quantitative anisotropy (QA) directly
    from the diffusion signal without assuming a parametric model. Native-space
    reconstruction is preferred for stroke populations with large lesions or ventricular
    enlargement: template-space reconstruction (QSDR) requires warping subject data into
    MNI space, and the registration degrades when the QA map deviates substantially from
    the healthy template -- as occurs with perilesional tissue loss and CSF expansion.
    Native-space GQI preserves the actual tissue signal and yields more reliable diffusion
    metric extraction. At AutoTrack time, DSI Studio performs its own on-the-fly
    registration of the atlas tracts into subject space. The FSL brain mask is applied
    during reconstruction if present.
    Skipped if sub-XXXX_eddy_unwarped_images.fib.gz already exists.

  Step 3 — AutoTrack:
    Runs atlas-based deterministic tractography on all default DSI Studio tract
    bundles using the ICBM152 template. Tract-specific parameters are set in the
    ADJUSTABLE PARAMETERS block at the top of the script.

  Step 4 — Statistics aggregation:
    Collects per-tract .stat.txt files from the AutoTrack output and combines
    them into a single TSV file with tracts listed alphabetically across columns
    and diffusion/morphological metrics as rows.

Outputs:
  sub-XXXX_eddy_unwarped_images.sz
  sub-XXXX_eddy_unwarped_images.fib.gz
  tract_stats/
    sub-XXXX_eddy_unwarped_images_autotrack.log
    sub-XXXX_eddy_unwarped_images_autotrack_stats.tsv
    sub-XXXX_eddy_unwarped_images_atk_work/

For session subjects, the session label is included in all output filenames:
  sub-XXXX_ses-XXXX_eddy_unwarped_images.sz  (etc.)

---------------------------------------
4. AutoTrack parameters
---------------------------------------

All parameters are set in the ADJUSTABLE PARAMETERS block at the top of
dwi_dsi_pipeline.sh. Defaults are listed below with brief rationale.

PARAM0 = 1.25
  GQI diffusion sampling length ratio. Standard in vivo default per DSI Studio docs.

TOLERANCE = 22,26,30
  Bundle recognition tolerance in mm. Comma-separated values trigger progressive
  retry per DSI Studio CLI behavior. For stroke populations with significant
  perilesional white matter distortion, this progressive scheme improves tract
  recovery on the lesioned side at the cost of slightly more false positives.

TRACK_VOXEL_RATIO = 2
  Controls the number of streamlines generated relative to tract volume.
  Higher values increase streamline counts but increase compute time.

CHECK_ENDING = 1
  Removes streamlines terminating in high-anisotropy regions (e.g., CSF/lesion
  borders). Enabled by default. DSI Studio automatically disables this for the
  cingulum regardless of this setting.

YIELD_RATE = 0.00001
  Early-termination threshold based on the ratio of accepted to attempted seeds.
  CLI default. Set to 0 to disable early stopping if low streamline counts are
  a concern.

TIP_ITERATION = 2
  Topology-informed pruning iterations. The DSI Studio developer documentation
  recommends 1-2 iterations. The GUI default is 4, but higher values increase
  false-negative rates in tracts with low streamline density. In stroke populations,
  perilesional tracts have systematically reduced streamline counts, making
  aggressive pruning likely to eliminate real surviving connections. Published
  DSI Studio AutoTrack studies in stroke populations use 2.
  Reference: Yeh et al. 2019, Neurotherapeutics (TIP original paper).

Thread count is not set explicitly; DSI Studio uses all available hardware threads,
or the value of NSLOTS if set (e.g., in a cluster job submission context).

Diffusion metrics embedded in the .fib.gz and available in tract statistics:
  fa    Fractional anisotropy
  md    Mean diffusivity
  ad    Axial diffusivity
  rd    Radial diffusivity
  rdi   Restricted diffusion imaging
  qa    Quantitative anisotropy (always computed by GQI regardless of other_output)

=======================================
III. BATCH PROCESSING
=======================================

Batch scripts run the single-subject scripts across all subjects in a cohort
parent directory. Both scripts search for dwi/ directories matching either the
single-timepoint or multi-timepoint structure under 2.preprocessing/.

---------------------------------------
1. Batch FSL preprocessing
---------------------------------------

Script: batch_dwi_fsl_preprocessing.sh

Usage:
  batch_dwi_fsl_preprocessing.sh /path/to/cohort_directory

Searches for:
  .../sub-XXXX/2.preprocessing/dwi/
  .../sub-XXXX/2.preprocessing/ses-XXXX/dwi/

Skips subjects where eddy_unwarped_images.nii.gz already exists.
Prints OK / SKIP / FAIL summary at completion.

---------------------------------------
2. Batch DSI Studio pipeline
---------------------------------------

Script: batch_dwi_dsi_processing.sh

Usage:
  batch_dwi_dsi_processing.sh /path/to/cohort_directory

Searches for:
  .../sub-XXXX/2.preprocessing/dwi/
  .../sub-XXXX/2.preprocessing/ses-XXXX/dwi/

Skips subjects where eddy_unwarped_images.nii.gz is absent (FSL not yet run).
Individual steps within the pipeline (SRC, FIB) are skipped per-subject if
their outputs already exist, allowing partial re-runs after interruption.
The AutoTrack workdir check will error if tract_stats already contains results
from a previous run -- move or remove the existing tract_stats directory first.

---------------------------------------
3. Expected folder structure
---------------------------------------

Single-timepoint subject:
  .../cohort/
    sub-XXXX/
      2.preprocessing/
        dwi/
          sub-XXXX_dir-AP_dwi.nii(.gz)
          sub-XXXX_dir-AP_dwi.bval
          ... (raw inputs)
          eddy_unwarped_images.nii.gz
          hifi_nodif_brain_mask.nii.gz
          ... (FSL outputs)
          sub-XXXX_eddy_unwarped_images.sz
          sub-XXXX_eddy_unwarped_images.fib.gz
          tract_stats/
            sub-XXXX_eddy_unwarped_images_autotrack_stats.tsv
            ...

Multi-timepoint subject:
  .../cohort/
    sub-XXXX/
      2.preprocessing/
        ses-XXXX/
          dwi/
            sub-XXXX_ses-XXXX_dir-AP_dwi.nii(.gz)
            sub-XXXX_ses-XXXX_dir-AP_dwi.bval
            ... (raw inputs)
            eddy_unwarped_images.nii.gz
            hifi_nodif_brain_mask.nii.gz
            ... (FSL outputs)
            sub-XXXX_ses-XXXX_eddy_unwarped_images.sz
            sub-XXXX_ses-XXXX_eddy_unwarped_images.fib.gz
            tract_stats/
              sub-XXXX_ses-XXXX_eddy_unwarped_images_autotrack_stats.tsv
              ...
        ses-YYYY/
          dwi/
            ...

=======================================
End of file
=======================================
