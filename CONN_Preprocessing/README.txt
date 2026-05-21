CONN Preprocessing Pipeline

This directory contains scripts used to generate CONN preprocessing configuration files (.cfg) and run CONN preprocessing for each subject.

The pipeline supports:
- PWA subjects (requires lesion mask)
- HC subjects (no lesion mask required, could also use on dementia subjects or others without lesion masking required)
- Single-subject processing
- Batch processing across a cohort
- Single or multiple sessions per subject

Scripts are located in:

SCC: /projectnb/openneuro-aphasia/ML_aphasia_scripts/CONN_Preprocessing/
Server: /Kiran/KiranLab2/ML-Aphasia project/4. Results, Analyses, Reports/ML_aphasia_scripts/CONN_Preprocessing/

--------------------------------------------------
DATA ORGANIZATION
--------------------------------------------------

Subjects should be organized as:

cohort_root/
  sub-XXXX/
    1.raw_data/
      (optional ses-YYYY/)
        DICOM/
          ...
    2.preprocessing/
      anat/
      func/
      lesion/      (PWA only)
      cfg/	(will be created by cfg builder script)

The 1.raw_data directory contains raw DICOM data downloaded from XNAT.

The 2.preprocessing directory contains converted NIfTI files and preprocessing outputs.


--------------------------------------------------
MULTI-SESSION SUBJECTS
--------------------------------------------------

If a subject has multiple sessions, session folders should be created inside 2.preprocessing.
The cfg/ directory is always placed directly inside 2.preprocessing, not inside individual
session folders:

sub-XXXX/
  2.preprocessing/
    ses-0001/
      anat/
      func/
      lesion/
    ses-0002/
      anat/
      func/
      lesion/
    cfg/

Session folders should follow the BIDS naming convention:

ses-0001
ses-0002

If no session folders exist, the scripts assume the subject has only one session.


--------------------------------------------------
EXPECTED FILE NAMING
--------------------------------------------------

Structural (required)

sub-XXXX[_ses-YYYY]_T1w.nii.gz

Functional (optional)

sub-XXXX[_ses-YYYY]_task-rest_bold.nii.gz

Lesion mask (PWA only)

sub-XXXX[_ses-YYYY]_desc-lesion_mask.nii.gz

Files must be placed inside:

anat/
func/
lesion/

directories within 2.preprocessing or 2.preprocessing/ses-*.


--------------------------------------------------
GENERATED CONFIGURATION FILES
--------------------------------------------------

Configuration files are written to:

2.preprocessing/cfg/

Examples:

sub-XXXX.cfg
sub-XXXX_ses-0001.cfg
sub-XXXX_ses-0002.cfg

Each .cfg file corresponds to one subject-session preprocessing run.


--------------------------------------------------
SCRIPTS OVERVIEW
--------------------------------------------------

CFG Builders (single subject) - Create CONN configuration files for a single subject.

build_conn_cfg_PWA_single.sh
build_conn_cfg_HC_single.sh


Batch CFG Builder - Runs the cfg builder across all subjects in a cohort.

batch_build_cfgs.sh


CONN Runner (single subject) - Runs CONN preprocessing using the cfg files generated for one subject.

run_conn_on_cfg.sh


Batch CONN Runner - Runs CONN preprocessing across all subjects in a cohort.

batch_run_conn_on_cfg.sh


--------------------------------------------------
SINGLE SUBJECT PROCESSING
--------------------------------------------------

Step 1 — Build configuration file


PWA subject: build_conn_cfg_PWA_single.sh /path/to/sub-XXXX


HC subject: build_conn_cfg_HC_single.sh /path/to/sub-XXXX


This script:

1. Detects whether 2.preprocessing exists.
2. Detects session folders (ses-*) if present.
3. Searches for required files in anat/, func/, and lesion/.
4. Generates .cfg files in:

2.preprocessing/cfg/


Example output:

sub-XXXX/2.preprocessing/cfg/sub-XXXX.cfg

or

sub-XXXX/2.preprocessing/cfg/sub-XXXX_ses-0001.cfg
sub-XXXX/2.preprocessing/cfg/sub-XXXX_ses-0002.cfg


ADJUSTABLE FIELDS

Each cfg builder script has an adjustable fields block near the top of the file:

DEFAULT_TR       - RepetitionTime in seconds, used if not found in the JSON sidecar.
                   Update this to match your acquisition before running.
DEFAULT_SLICEORDER - Slice acquisition order written to the .cfg file.
                   Update this to match your acquisition (e.g. "interleaved (Siemens)",
                   "ascending", "descending"). All files will be written using the slice order specified here, so if running batch with multiple orders, manually update .cfg files after building


--------------------------------------------------
STEP 2 — RUN CONN PREPROCESSING
--------------------------------------------------

run_conn_on_cfg.sh /path/to/sub-XXXX


This script:

1. Locates the subject's cfg directory.
2. Finds all .cfg files.
3. Runs CONN preprocessing on each file.

Example:

sub-XXXX/2.preprocessing/cfg/sub-XXXX.cfg

or multiple sessions:

sub-XXXX_ses-0001.cfg
sub-XXXX_ses-0002.cfg

Each configuration is processed sequentially.

NOTE: The MATLAB, SPM, and CONN module paths and versions in run_conn_on_cfg.sh are
specific to the BU Shared Computing Cluster (SCC). Users on other systems must update
these paths to match their local installations before running.


--------------------------------------------------
BATCH PROCESSING
--------------------------------------------------

Batch scripts allow processing entire cohorts.


--------------------------------------------------
BATCH STEP 1 — BUILD CFG FILES
--------------------------------------------------

PWA cohort

batch_build_cfgs.sh PWA /path/to/cohort_root


HC cohort

batch_build_cfgs.sh HC /path/to/cohort_root


This script:

1. Searches for directories matching

sub-*

2. Calls the appropriate single-subject cfg builder for each subject.


--------------------------------------------------
BATCH STEP 2 — RUN CONN PREPROCESSING
--------------------------------------------------

batch_run_conn_on_cfg.sh /path/to/cohort_root


This script:

1. Finds every sub-* directory.
2. Calls run_conn_on_cfg.sh for each subject.
3. Runs preprocessing for every .cfg file found.
4. Skips subjects where all expected output files already exist (resumable).


--------------------------------------------------
LESION WARPING
--------------------------------------------------

warp_lesion_to_target.sh warps a lesion mask from one structural image space to another
using ANTs SyN registration. This is used when a lesion mask exists for one session and
needs to be propagated to a second session.

Usage:

bash warp_lesion_to_target.sh <source_lesion> <source_anat> <target_anat>

Inputs may be .nii or .nii.gz. Output is always written as uncompressed .nii.

The source anat is registered to the target anat using cost-function masking: the lesion
is passed as a moving-image mask to exclude infarcted voxels from the registration metric.
The resulting deformation field is then applied to the lesion mask to warp it into the
target space.

Output is written to the lesion/ directory adjacent to the target anat/:

sub-XXXX/2.preprocessing/ses-YYYY/lesion/sub-XXXX_ses-YYYY_desc-lesion_mask.nii

For batch processing across a cohort:

bash batch_warp_lesion.sh <parent_dir> <input_session> <target_session>

NOTE: warp_lesion_to_target.sh requires ANTs and FSL. Module load commands are specific
to the BU SCC and must be updated for other systems.


--------------------------------------------------
EXAMPLE WORKFLOW
--------------------------------------------------

Single subject

build_conn_cfg_PWA_single.sh sub-BUBA001
run_conn_on_cfg.sh sub-BUBA001


Entire cohort

batch_build_cfgs.sh PWA /projectnb/openneuro-aphasia/CBR_Cohort_B
batch_run_conn_on_cfg.sh /projectnb/openneuro-aphasia/CBR_Cohort_B


--------------------------------------------------
NOTES
--------------------------------------------------

Functional scans are optional. If absent, structural preprocessing will run only.

Lesion masks are required for PWA subjects.

All preprocessing outputs are written by CONN into the corresponding anat/, lesion/, and
func/ directories.

DEFAULT_TR and DEFAULT_SLICEORDER in the cfg builder scripts must be verified and updated
to match your acquisition before running. TR is read automatically from JSON sidecars when
present; slice order is always taken from the DEFAULT_SLICEORDER field.

Module paths in run_conn_on_cfg.sh and warp_lesion_to_target.sh are specific to the BU SCC.
Update these before use on other systems.


--------------------------------------------------
SUMMARY
--------------------------------------------------

This pipeline allows flexible preprocessing across cohorts with:

- optional sessions
- optional functional scans
- lesion-aware normalization for PWA
- cost-function masked structural registration for lesion warping
- standardized organization of cfg files
- support for both single-subject and batch processing