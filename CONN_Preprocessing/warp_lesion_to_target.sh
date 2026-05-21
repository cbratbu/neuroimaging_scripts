#!/bin/bash
# Warp a lesion mask from one structural image space to a target structural image space using ANTs SyN
# Usage: bash warp_lesion_to_target.sh <source_lesion> <source_anat> <target_anat>
# Inputs may be .nii or .nii.gz
# Sub/session labels parsed from BIDS-convention filenames (sub-[ID]_ses-[session]_...)
# Output: sub-[ID]_ses-[session]_desc-lesion_mask.nii saved to lesion/ adjacent to target anat/
# Requires ANTs and FSL modules

# ----- Inputs -----
SOURCE_LESION="$1"
SOURCE_ANAT="$2"
TARGET_ANAT="$3"

ANTS_THREADS=$(nproc)

# ----- Usage check -----
if [ -z "$SOURCE_LESION" ] || [ -z "$SOURCE_ANAT" ] || [ -z "$TARGET_ANAT" ]; then
    echo "USAGE ERROR: bash warp_lesion_to_target.sh <source_lesion> <source_anat> <target_anat>"
    exit 1
fi

# ----- File existence check -----
for f in "$SOURCE_LESION" "$SOURCE_ANAT" "$TARGET_ANAT"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

# ----- Parse sub/session from target anat filename -----
TARGET_BASENAME=$(basename "$TARGET_ANAT")
TARGET_BASENAME="${TARGET_BASENAME%.nii.gz}"
TARGET_BASENAME="${TARGET_BASENAME%.nii}"

TARGET_SUB=$(echo "$TARGET_BASENAME" | grep -oP 'sub-[^_]+')
TARGET_SES=$(echo "$TARGET_BASENAME" | grep -oP 'ses-[^_]+')

if [ -z "$TARGET_SUB" ]; then
    echo "ERROR: Could not parse subject label from target filename: $(basename $TARGET_ANAT)"
    exit 1
fi
if [ -z "$TARGET_SES" ]; then
    echo "ERROR: Could not parse session label from target filename: $(basename $TARGET_ANAT)"
    echo "ERROR: Session label (ses-[session]) is required in all filenames"
    exit 1
fi

# ----- Derive output path from target anat location -----
TARGET_ANAT_DIR=$(dirname "$TARGET_ANAT")
TARGET_SES_DIR=$(dirname "$TARGET_ANAT_DIR")

OUTPUT_DIR="${TARGET_SES_DIR}/lesion"
OUTPUT_NAME="${TARGET_SUB}_${TARGET_SES}_desc-lesion_mask.nii"
OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_NAME}"

mkdir -p "$OUTPUT_DIR"

# ----- Temp dir for ANTs intermediate files -----
TMP_DIR=$(mktemp -d)
WARP_PREFIX="${TMP_DIR}/ants_"

# ----- Load modules -----
module load ants/2.6.2
module load fsl

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$ANTS_THREADS

# ----- Register source anat to target anat -----
# lesion passed as moving-image mask to exclude infarcted voxels
# from the registration metric, preventing the deformation field from being driven by
# non-corresponding tissue. The resulting warp is then applied to the lesion separately.
antsRegistrationSyN.sh \
    -d 3 \
    -f "$TARGET_ANAT" \
    -m "$SOURCE_ANAT" \
    -x "NULL,$SOURCE_LESION" \
    -o "$WARP_PREFIX" \
    -t s \
    -n "$ANTS_THREADS"

if [ $? -ne 0 ]; then
    echo "ERROR: ANTs registration failed for ${TARGET_SUB} ${TARGET_SES}"
    rm -rf "$TMP_DIR"
    exit 1
fi

# ----- Apply warp to lesion mask -----
antsApplyTransforms \
    -d 3 \
    -i "$SOURCE_LESION" \
    -r "$TARGET_ANAT" \
    -o "${TMP_DIR}/lesion_warped.nii" \
    -t "${WARP_PREFIX}1Warp.nii.gz" \
    -t "${WARP_PREFIX}0GenericAffine.mat" \
    -n NearestNeighbor

if [ $? -ne 0 ]; then
    echo "ERROR: antsApplyTransforms failed for ${TARGET_SUB} ${TARGET_SES}"
    rm -rf "$TMP_DIR"
    exit 1
fi

# ----- Binarize at midpoint of voxel intensity range -----
# is used as a safeguard for any non-binary input masks
VOXEL_MIN=$(fslstats "${TMP_DIR}/lesion_warped.nii" -R | awk '{print $1}')
VOXEL_MAX=$(fslstats "${TMP_DIR}/lesion_warped.nii" -R | awk '{print $2}')
THRESHOLD=$(echo "($VOXEL_MIN + $VOXEL_MAX) / 2" | bc -l)

FSLOUTPUTTYPE=NIFTI fslmaths "${TMP_DIR}/lesion_warped.nii" -thr "$THRESHOLD" -bin "$OUTPUT_PATH"

if [ $? -ne 0 ]; then
    echo "ERROR: Binarization failed for ${TARGET_SUB} ${TARGET_SES}"
    rm -rf "$TMP_DIR"
    exit 1
fi

# ----- Cleanup -----
rm -rf "$TMP_DIR"

echo "SUCCESS: ${TARGET_SUB} ${TARGET_SES} -> ${OUTPUT_PATH}"
exit 0