#!/bin/bash

# -----------------------------
# make_acqparams_index.sh
# -----------------------------
# Usage:
#   ./make_acqparams_index.sh /path/to/dwi_folder
#
# Expects inside folder:
#   * *_dir-AP_dwi.json
#   * *_dir-PA_dwi.json
#   * *_dir-AP_dwi.nii (for volume count)
# -----------------------------

set -e  # exit if any command fails

if [ -z "$1" ]; then
    echo "ERROR: No directory supplied."
    echo "Usage: $0 /path/to/dwi_folder"
    exit 1
fi

DWI_DIR="$1"
cd "$DWI_DIR" || { echo "Directory not found."; exit 1; }

AP_JSON=$(ls *dir-AP_dwi.json 2>/dev/null || true)
PA_JSON=$(ls *dir-PA_dwi.json 2>/dev/null || true)
AP_NII=$(ls *dir-AP_dwi.nii* | head -n 1)

if [ ! -f "$AP_JSON" ] || [ ! -f "$PA_JSON" ]; then
    echo "ERROR: AP or PA JSON file missing."
    exit 1
fi

echo "Using:"
echo "  AP JSON: $AP_JSON"
echo "  PA JSON: $PA_JSON"
echo "  AP NIfTI for index count: $AP_NII"
echo ""

# -----------------------
# Get values from JSONs
# -----------------------
get_pe_vec () {
    case "$1" in
        i) echo "1 0 0" ;;
        i-) echo "-1 0 0" ;;
        j) echo "0 1 0" ;;
        j-) echo "0 -1 0" ;;
        k) echo "0 0 1" ;;
        k-) echo "0 0 -1" ;;
        *) echo "Unknown PhaseEncodingDirection: $1" && exit 1 ;;
    esac
}

AP_PE=$(grep -o '"PhaseEncodingDirection": *"[^"]*"' "$AP_JSON" | awk -F'"' '{print $4}')
AP_TRO=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$AP_JSON" | awk '{print $2}')

PA_PE=$(grep -o '"PhaseEncodingDirection": *"[^"]*"' "$PA_JSON" | awk -F'"' '{print $4}')
PA_TRO=$(grep -o '"TotalReadoutTime": *[0-9.]*' "$PA_JSON" | awk '{print $2}')

AP_VEC=$(get_pe_vec "$AP_PE")
PA_VEC=$(get_pe_vec "$PA_PE")

# -----------------------
# Write acqparams.txt
# -----------------------
echo "Writing acqparams.txt ..."
echo "$AP_VEC $AP_TRO" > acqparams.txt
echo "$PA_VEC $PA_TRO" >> acqparams.txt

echo "acqparams.txt created:"
cat acqparams.txt
echo ""

# -----------------------
# Create index.txt
# -----------------------
echo "Creating index.txt ..."
NVOLS=$(fslval "$AP_NII" dim4)

echo -n "" > index.txt
for ((i=1; i<=NVOLS; i++)); do
    echo -n "1 " >> index.txt
done
echo "" >> index.txt

echo "index.txt created (showing first 10 values):"
cat index.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10"..."}'
echo ""

echo "Done!"

