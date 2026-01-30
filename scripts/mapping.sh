#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
    echo "Usage: $0 -r REF_DIR -f REF_NAME -o OUT_DIR -q FASTQ_DIR -n VIRUS_NAME -i SAMPLE_PREFIX"
}

#######################
# Parse arguments
#######################
while getopts "r:f:o:q:n:i:h" opt; do
    case $opt in
        r) REF_DIR="$OPTARG" ;;
        f) REF_NAME="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        q) FASTQ_DIR="$OPTARG" ;;
        n) VIRUS_NAME="$OPTARG" ;;
        i) SAMPLE_PREFIX="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check required arguments
if [[ -z "${REF_DIR:-}" || -z "${REF_NAME:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_DIR:-}" || -z "${VIRUS_NAME:-}" || -z "${SAMPLE_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# Step 1: Add LTR to reference
#######################
echo "Add LTR sequences to the host reference genome"

MASKED_LTR5="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR5_withprimer.fa"
MASKED_LTR3="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR3_withprimer.fa"

if [[ ! -f "$MASKED_LTR3" ]]; then
    echo "Creating $MASKED_LTR3"
    cat "${REF_DIR}/$REF_NAME" \
        "${REF_DIR}/${VIRUS_NAME}_endU3RU5_withprimer.fa" \
        > "$MASKED_LTR3"
else
    echo "Found existing $MASKED_LTR3 — skipping"
fi

if [[ ! -f "$MASKED_LTR5" ]]; then
    echo "Creating $MASKED_LTR5"
    cat "${REF_DIR}/$REF_NAME" \
        "${REF_DIR}/${VIRUS_NAME}_startU3_withprimer.fa" \
        > "$MASKED_LTR5"
else
    echo "Found existing $MASKED_LTR5 — skipping"
fi


#######################
# Step 2: Index hybrid references
#######################
MMI_LTR3="${MASKED_LTR3%.fa}.mmi"
MMI_LTR5="${MASKED_LTR5%.fa}.mmi"

if [[ ! -f "$MMI_LTR3" ]]; then
    echo "Indexing $MASKED_LTR3"
    minimap2 -d "$MMI_LTR3" "$MASKED_LTR3"
else
    echo "Found existing $MMI_LTR3 — skipping"
fi

if [[ ! -f "$MMI_LTR5" ]]; then
    echo "Indexing $MASKED_LTR5"
    minimap2 -d "$MMI_LTR5" "$MASKED_LTR5"
else
    echo "Found existing $MMI_LTR5 — skipping"
fi

#######################
# Step 3: Map reads for each LTR
#######################
for LTR_NUM in 3 5; do
    if [ "$LTR_NUM" -eq 3 ]; then
        MMI_FILE="${MASKED_LTR3%.fa}.mmi"
    else
        MMI_FILE="${MASKED_LTR5%.fa}.mmi"
    fi

    FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_filtered_size_SUP.fastq"
    SAM_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.sam"
    
    echo "Mapping $FASTQ_FILE to LTR${LTR_NUM} reference..."
    minimap2 -ax map-ont -t 8 "$MMI_FILE" "$FASTQ_FILE" > "$SAM_FILE"

    # Convert to PAF
    PAF_DIR="${OUT_DIR}/paf"
    mkdir -p "$PAF_DIR"
    paftools.js sam2paf "$SAM_FILE" > "${PAF_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.paf"

done

#######################
# Step 4: Cleanup intermediate files
#######################
echo "####### Cleaning up intermediate files..."
find "$OUT_DIR" -type f ! -name "*.paf" -delete

echo "####### DONE - $SAMPLE_PREFIX"
