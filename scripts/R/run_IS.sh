#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
    echo "Usage: $0 -i SAMPLE_NAME -r R_PACKAGE_PATH -o OUT_PATH -p INPUT_PAF_PATH -u INPUT_UMI_PATH -a ASSEMBLY \\"
    echo "          -5 TNAME_LTR5 -3 TNAME_LTR3 \\"
    echo "          -g MAXGAP_IS -q MAPQ -n NB_READ_PER_UMI -s MAXGAP_SHS -m MMS -t THRESHOLD_RAW -w WIN_MERGE"
}

#######################
# Parse arguments
#######################
while getopts "i:r:o:p:u:a:5:3:g:q:n:s:m:t:w:h" opt; do
    case $opt in
        i) SAMPLE_NAME="$OPTARG" ;;
        r) R_PACKAGE_PATH="$OPTARG" ;;
        o) OUT_PATH="$OPTARG" ;;
        p) INPUT_PAF_PATH="$OPTARG" ;;
        u) INPUT_UMI_PATH="$OPTARG" ;;
        a) ASSEMBLY="$OPTARG" ;;
        5) TARGET_LTR5="$OPTARG" ;;
        3) TARGET_LTR3="$OPTARG" ;;
        g) MAXGAP_IS="$OPTARG" ;;
        q) MAPQ_VAL="$OPTARG" ;;
        n) NB_READ_PER_UMI="$OPTARG" ;;
        s) MAXGAP_SHS="$OPTARG" ;;
        m) MMS="$OPTARG" ;;
        t) THRESHOLD_RAW="$OPTARG" ;;
        w) WIN_MERGE="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

#######################
# Check required arguments
#######################
if [[ -z "${SAMPLE_NAME:-}" || -z "${R_PACKAGE_PATH:-}" || -z "${OUT_PATH:-}" || -z "${INPUT_PAF_PATH:-}" || -z "${INPUT_UMI_PATH:-}" || -z "${ASSEMBLY:-}" || \
      -z "${TARGET_LTR5:-}" || -z "${TARGET_LTR3:-}" || \
      -z "${MAXGAP_IS:-}" || -z "${MAPQ_VAL:-}" || -z "${NB_READ_PER_UMI:-}" || -z "${MAXGAP_SHS:-}" || \
      -z "${MMS:-}" || -z "${THRESHOLD_RAW:-}" || -z "${WIN_MERGE:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

#######################
# Pre-run cleanup: remove existing files with same sample prefix
#######################
if [[ -d "$OUT_PATH" ]]; then
    echo "Checking for existing files with prefix ${SAMPLE_NAME} in ${OUT_PATH}..."

    shopt -s nullglob
    existing_files=("${OUT_PATH}${SAMPLE_NAME}"*)
    shopt -u nullglob

    if (( ${#existing_files[@]} > 0 )); then
        echo "Found existing files. Removing them..."
        rm -f "${existing_files[@]}"
        echo "Previous files removed."
    else
        echo "No existing files found."
    fi
fi

mkdir -p "$OUT_PATH"

R_SCRIPT_DIR="$(cd "$(dirname "$R_PACKAGE_PATH")" && pwd)"

#######################
# Step 1: First R script
#######################
echo "Starting first R script..."

Rscript "${R_SCRIPT_DIR}/run_IS_script_1.R" \
    "$R_PACKAGE_PATH" \
    "$SAMPLE_NAME" \
    "$OUT_PATH" \
    "$INPUT_PAF_PATH" \
    "$INPUT_UMI_PATH" \
    "$TARGET_LTR5" \
    "$TARGET_LTR3" \
    "$MAPQ_VAL" \
    "$ASSEMBLY"

echo "First R script completed."

#######################
# Step 2: UMI clustering
#######################
echo "Starting UMI clustering..."

for file in ${OUT_PATH}${SAMPLE_NAME}*data2*LTR*.txt; do
    python3 ${R_SCRIPT_DIR}/UMI_clustering_hamming_ref.py \
        "$file" \
        "${file}_hamming_ref_mms${MMS}.txt" \
        --mismatch_threshold "$MMS"
done

echo "UMI clustering completed."

#######################
# Step 3: Merge chromosome files
#######################
echo "Merging LTR5 files..."
awk 'FNR==1 && NR!=1 { next } { print }' \
    ${OUT_PATH}${SAMPLE_NAME}*data2*_LTR5.txt_hamming_ref_mms${MMS}.txt \
    > "${OUT_PATH}${SAMPLE_NAME}_merged_LTR5_mms${MMS}.txt"

echo "Merging LTR3 files..."
awk 'FNR==1 && NR!=1 { next } { print }' \
    ${OUT_PATH}${SAMPLE_NAME}*data2*_LTR3.txt_hamming_ref_mms${MMS}.txt \
    > "${OUT_PATH}${SAMPLE_NAME}_merged_LTR3_mms${MMS}.txt"

#######################
# Step 4: Second R script
#######################
echo "Starting second R script..."

Rscript "${R_SCRIPT_DIR}/run_IS_script_2.R" \
    "$R_PACKAGE_PATH" \
    "$SAMPLE_NAME" \
    "$OUT_PATH" \
    "$MAXGAP_IS" \
    "$NB_READ_PER_UMI" \
    "$MAXGAP_SHS" \
    "$MMS" \
    "$WIN_MERGE" \
    "$THRESHOLD_RAW"

#######################
# Step 5: Cleanup intermediate files
#######################
echo "####### Cleaning up intermediate files..."
find "$OUT_PATH" -type f ! -name "*clonality*" -delete

echo "####### JOB FINISHED SUCCESSFULLY - $SAMPLE_NAME"
