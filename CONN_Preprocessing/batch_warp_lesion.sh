#!/bin/bash
# Batch wrapper to warp lesion masks from an input session to a target session across subjects
# Usage: bash batch_warp_lesion.sh <parent_dir> <input_session> <target_session>
# Expects standard 2.preprocessing directory structure and T1w naming conventions per subject
# Inputs may be .nii or .nii.gz; output is always written as .nii by warp_lesion_to_target.sh
# Calls: /projectnb/openneuro-aphasia/ML_aphasia_scripts/CONN_Preprocessing/warp_lesion_to_target.sh

WARP_SCRIPT="/projectnb/openneuro-aphasia/ML_aphasia_scripts/CONN_Preprocessing/warp_lesion_to_target.sh"

PARENT_DIR="$1"
INPUT_SES="$2"
TARGET_SES="$3"

if [ -z "$PARENT_DIR" ] || [ -z "$INPUT_SES" ] || [ -z "$TARGET_SES" ]; then
    echo "Usage: bash batch_warp_lesion.sh <parent_dir> <input_session> <target_session>"
    exit 1
fi

if [ ! -d "$PARENT_DIR" ]; then
    echo "ERROR: Parent directory not found: $PARENT_DIR"
    exit 1
fi

if [ ! -f "$WARP_SCRIPT" ]; then
    echo "ERROR: Warp script not found: $WARP_SCRIPT"
    exit 1
fi

# Helper: return path if .nii or .nii.gz exists, otherwise empty string
find_nii() {
    local base="$1"
    if [ -f "${base}.nii" ]; then
        echo "${base}.nii"
    elif [ -f "${base}.nii.gz" ]; then
        echo "${base}.nii.gz"
    else
        echo ""
    fi
}

# ----- Counters -----
n_success=0
n_fail=0
n_skip_missing=0
n_skip_complete=0

declare -a skipped_missing
declare -a skipped_complete
declare -a failed

# ----- Loop over subjects -----
for sub_dir in "$PARENT_DIR"/sub-*/; do
    sub_id=$(basename "$sub_dir")

    source_lesion=$(find_nii "${sub_dir}2.preprocessing/${INPUT_SES}/lesion/${sub_id}_${INPUT_SES}_desc-lesion_mask")
    source_anat=$(find_nii "${sub_dir}2.preprocessing/${INPUT_SES}/anat/${sub_id}_${INPUT_SES}_T1w")
    target_anat=$(find_nii "${sub_dir}2.preprocessing/${TARGET_SES}/anat/${sub_id}_${TARGET_SES}_T1w")

    # Output is always .nii (warp_lesion_to_target.sh writes uncompressed)
    expected_output="${sub_dir}2.preprocessing/${TARGET_SES}/lesion/${sub_id}_${TARGET_SES}_desc-lesion_mask.nii"

    # ----- Skip if already complete -----
    if [ -f "$expected_output" ]; then
        echo "ALREADY COMPLETE (lesion mask exists in target session): $sub_id"
        n_skip_complete=$((n_skip_complete + 1))
        skipped_complete+=("$sub_id")
        continue
    fi

    # ----- Skip if any required file is missing -----
    missing=()
    [ -z "$source_lesion" ] && missing+=("source lesion: ${sub_dir}2.preprocessing/${INPUT_SES}/lesion/${sub_id}_${INPUT_SES}_desc-lesion_mask.nii[.gz]")
    [ -z "$source_anat" ]   && missing+=("source anat: ${sub_dir}2.preprocessing/${INPUT_SES}/anat/${sub_id}_${INPUT_SES}_T1w.nii[.gz]")
    [ -z "$target_anat" ]   && missing+=("target anat: ${sub_dir}2.preprocessing/${TARGET_SES}/anat/${sub_id}_${TARGET_SES}_T1w.nii[.gz]")

    if [ ${#missing[@]} -gt 0 ]; then
        echo "SKIP (missing files): $sub_id"
        for m in "${missing[@]}"; do
            echo "  - $m"
        done
        n_skip_missing=$((n_skip_missing + 1))
        skipped_missing+=("$sub_id")
        continue
    fi

    # ----- Call single-subject warp script -----
    bash "$WARP_SCRIPT" "$source_lesion" "$source_anat" "$target_anat"

    if [ $? -eq 0 ]; then
        n_success=$((n_success + 1))
    else
        echo "FAILED: $sub_id"
        n_fail=$((n_fail + 1))
        failed+=("$sub_id")
    fi

done

# ----- Summary -----
echo ""
echo "===== Batch Warp Summary ====="
echo "Input session:  $INPUT_SES"
echo "Target session: $TARGET_SES"
echo ""
echo "Successful:              $n_success"
echo "Failed:                  $n_fail"
echo "Skipped (missing files): $n_skip_missing"
echo "Skipped (already done):  $n_skip_complete"

if [ ${#failed[@]} -gt 0 ]; then
    echo ""
    echo "Failed subjects:"
    for s in "${failed[@]}"; do echo "  $s"; done
fi

if [ ${#skipped_missing[@]} -gt 0 ]; then
    echo ""
    echo "Skipped (missing files):"
    for s in "${skipped_missing[@]}"; do echo "  $s"; done
fi

if [ ${#skipped_complete[@]} -gt 0 ]; then
    echo ""
    echo "Skipped (already complete):"
    for s in "${skipped_complete[@]}"; do echo "  $s"; done
fi

exit 0