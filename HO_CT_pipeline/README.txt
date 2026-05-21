README — HO_CT_pipeline
========================

Surface-based cortical thickness extraction in Harvard-Oxford ROIs using FreeSurfer recon-all outputs. Per-ROI mean thickness is extracted via surface-based atlas transfer (mri_label2label + mris_anatomical_stats), avoiding volumetric registration to MNI space.


Dependencies
------------
FreeSurfer 6.0, FSL (setup step only), Python 3
Optional: scipy (for Shapiro-Wilk normality test in normative aggregation; falls back to NaN if unavailable)
GNU parallel (for HCP normative batch processing only)


Atlas
-----
Harvard-Oxford cortical atlas, 48 ROIs, thr25 maxprob (FSL distribution).
Makris et al., Schizophr Res. 2006;83(2-3):155-171.
Desikan et al., NeuroImage. 2006;31(3):968-980.

The atlas NIfTI files are distributed with FSL and are not included in this
repository. They are read directly from $FSLDIR/data/atlases/HarvardOxford/
during the one-time setup step.


Reference files (included)
--------------------------
HO48_1based.lut               FreeSurfer color table for HO atlas (48 ROIs, 1-based indexing)
HO48_labels_1based.tsv        Index-to-name mapping (tab-separated, no header, col 0 = index, col 1 = display name)

These files define the ROI identity mapping used throughout the pipeline. Do not
modify them unless you carefully verify that all index-to-name mappings remain
consistent across the LUT, TSV, and any normative reference data.


Normative data (not included)
------------------------------
normative_thickness_ROI_stats.csv    Per-ROI normative statistics (mean, SD, median, IQR) derived from HCP S1200
normative_subject_volume_means.csv   Normative brain volume reference

These files are not included in this repository because they are derived from
the HCP S1200 dataset, which requires data use through HCP
(https://db.humanconnectome.org). Scripts to regenerate them from a local HCP
S1200 download are provided — see the Normative Data Generation section below.
Patient-vs-norms comparisons in cortical_thickness_metrics_pipeline.sh will be
skipped if these files are absent.


Setup (one-time, required before any subject processing)
---------------------------------------------------------
Projects the HO atlas from MNI152 volumetric space onto the fsaverage surface
and builds annotation files. Must be run once before processing any subjects.
Outputs are reused by all subsequent subject-level runs.

  module load freesurfer/6.0 fsl
  bash ho_ct_setup_fsaverage.sh

This script:
  - Copies fsaverage to a writable location (atlas/subjects_tmp/fsaverage/)
    because mris_label2annot writes output to SUBJECTS_DIR and the FreeSurfer
    installation directory is read-only on the SCC
  - Projects both thr0 and thr25 HO atlas variants onto the fsaverage surface
    using mri_vol2surf with --mni152reg and nearest-neighbor interpolation
  - Extracts per-label surface label files for each of 48 ROIs per hemisphere
  - Assembles per-hemisphere .annot files using mris_label2annot

Outputs saved to atlas/:
  lh.ho_thr25.annot, rh.ho_thr25.annot     fsaverage HO annotations (thr25, recommended)
  lh.ho_thr0.annot,  rh.ho_thr0.annot      fsaverage HO annotations (thr0, available if needed)
  fsaverage_labels/thr25/                   Per-label source files used by per-subject processing
  fsaverage_labels/thr0/                    Per-label source files (thr0 variant)
  subjects_tmp/fsaverage/                   Writable fsaverage copy


Per-subject processing
-----------------------
  module load freesurfer/6.0
  bash cortical_thickness_metrics_pipeline.sh /path/to/<freesurfer_subject_dir> \
    out=/path/to/output_dir [thr=thr25|thr0]

The argument must be the FreeSurfer subject directory containing mri/, surf/,
label/, stats/ (i.e., the recon-all output root for that subject, typically
the directory named sub-XXXX produced by recon-all). Default atlas threshold
is thr25.

This script:
  Step 1 — Transfers HO atlas labels from fsaverage to the subject's native
            surface via mri_label2label using sulcal topology alignment
            (--regmethod surface). No volumetric MNI registration is performed.
  Step 2 — Assembles a subject-space HO annotation file with mris_label2annot.
  Step 3 — Extracts per-ROI mean cortical thickness, SD, vertex count, and gray
            matter volume with mris_anatomical_stats, restricted to cortical
            vertices via the subject's cortex.label.
  Step 4 — Parses stats and writes thickness CSV with bilateral mean (Mean_CT),
            laterality index (LI), and within-subject z-scores.
  Step 5 — Extracts eTIV, total gray matter, white matter, and CSF volumes from
            aseg.stats.
  Step 6 — Compares subject volumes to normative reference (skipped if
            normative_subject_volume_means.csv is absent).
  Step 7 — Compares subject cortical thickness to normative reference (skipped
            if normative_thickness_ROI_stats.csv is absent).

All intermediate files are written to a per-run working directory
(<out_prefix>_cortical_thickness_pipeline/work/) so the subject's recon-all
outputs are not modified.


Outputs (written to out= directory)
-------------------------------------
All outputs are placed inside a subdirectory named <prefix>_cortical_thickness_pipeline/:

  <prefix>_HO48_thickness.csv
    Per-ROI mean cortical thickness (LH, RH), SD, vertex count, gray matter
    volume, bilateral mean (Mean_CT), laterality index (LI), and within-subject
    z-scores relative to global and hemispheric ROI distributions.

  <prefix>_subject_volumes.csv
    eTIV, total gray matter, cerebral white matter, and CSF volumes from aseg.stats.

  <prefix>_patient_vs_norms.csv
    Per-ROI z-scores and robust z-scores (IQR-based) relative to normative
    healthy control data. Requires normative_thickness_ROI_stats.csv.

  <prefix>_volumes_vs_norms.csv
    Volume z-scores relative to normative reference. Requires
    normative_subject_volume_means.csv.

  work/
    Working directory containing intermediate files, transferred labels,
    subject-space annotation, and stats files.


Normative data generation (HCP S1200)
---------------------------------------
The normative reference CSVs are generated by processing all valid HCP S1200
subjects through the same surface-based pipeline and aggregating the results.
This requires a local HCP S1200 download with data use agreement
(https://db.humanconnectome.org). Three scripts handle this:

  hcp_normative_ct_subject.sh
    Processes a single HCP subject. Expects the HCP directory structure:
      <hcp_base>/<subject_id>/T1w/<subject_id>/surf/
    Performs label transfer, annotation assembly, mris_anatomical_stats, and
    volume extraction, writing per-subject CSVs to <output_base>/<subject_id>/.
    Skips subjects already complete (both _HO48_thickness.csv and
    _subject_volumes.csv present). Can be run standalone or called by the batch
    script. Hardcoded to thr25.

    Usage:
      module load freesurfer/6.0
      bash hcp_normative_ct_subject.sh <subject_id> <hcp_base_dir> <output_base_dir>

  hcp_normative_ct_batch.sh
    Runs hcp_normative_ct_subject.sh across all valid HCP subjects in parallel
    using GNU parallel, then calls the aggregation script. Adjustable fields at
    the top of the script control paths and job count.

    Usage:
      module load freesurfer/6.0 parallel
      bash hcp_normative_ct_batch.sh [--jobs N] [--aggregate-only]

    Options:
      --jobs N           Number of parallel jobs (default: 16)
      --aggregate-only   Skip per-subject processing; run aggregation only on
                         existing CSVs (useful after interrupted runs)

    The batch script uses GNU parallel with --resume-failed and a joblog, so
    interrupted runs can be safely restarted. Already-completed subjects are
    skipped both by the joblog and by the per-subject skip check.

  hcp_normative_ct_aggregate.py
    Aggregates per-subject _HO48_thickness.csv files into a single normative
    statistics CSV. Computes per-ROI mean, SD, median, IQR, and Shapiro-Wilk
    normality p-value across all subjects for LH, RH, bilateral mean (Both),
    and laterality index (LI). Requires scipy for Shapiro-Wilk; falls back to
    NaN if unavailable.

    Usage:
      python3 hcp_normative_ct_aggregate.py <output_base_dir> <norms_out_csv>

    <output_base_dir> should be the same directory passed to the batch script,
    containing one subdirectory per subject. The script discovers CSVs via
    glob pattern */*_HO48_thickness.csv.

Adjustable fields in hcp_normative_ct_batch.sh (update before running):
  HCP_BASE          Path to local HCP S1200 download
  OUT_BASE          Output directory for per-subject CSVs
  PIPELINE_DIR      Path to this pipeline directory (for writing final norms CSV)
  N_JOBS            Number of parallel workers (default: 16)

The final normative CSV is written to:
  <PIPELINE_DIR>/normative_thickness_ROI_stats.csv

This file must be regenerated if the atlas threshold or pipeline method changes,
as it must be produced by the same surface-based procedure applied to subjects.


Notes
------
- All SCC module load commands (freesurfer/6.0, fsl, parallel) are specific to
  the BU Shared Computing Cluster and must be updated for other systems.
- The HCP FreeSurfer outputs used here were produced by the HCP minimal
  preprocessing pipeline (Glasser et al., 2013), which uses a different
  surface registration target than standard recon-all. Label transfer via
  sphere.reg is valid in both cases as both register to the fsaverage sphere.
- thr25 is recommended over thr0 for most analyses; thr0 includes all atlas
  voxels regardless of overlap probability, increasing the chance of ambiguous
  border assignments.


References
----------
Fischl B. FreeSurfer. NeuroImage. 2012;62(2):774-781.
  doi:10.1016/j.neuroimage.2012.01.021

Fischl B, Dale AM. Measuring the thickness of the human cerebral cortex from
  magnetic resonance images. PNAS. 2000;97(20):11050-11055.

Glasser MF et al. The minimal preprocessing pipelines for the Human Connectome
  Project. NeuroImage. 2013;80:105-124. doi:10.1016/j.neuroimage.2013.04.127

Makris N et al. Decreased volume of left and total anterior insular lobule in
  schizophrenia. Schizophr Res. 2006;83(2-3):155-171.

Desikan RS et al. An automated labeling system for subdividing the human cerebral
  cortex on MRI scans into gyral based regions of interest. NeuroImage.
  2006;31(3):968-980.