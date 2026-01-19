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

#######################
# Check required arguments
#######################
if [[ -z "${REF_DIR:-}" || -z "${REF_NAME:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_DIR:-}" || -z "${VIRUS_NAME:-}" || -z "${SAMPLE_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# Step 1: Add LTR to reference
#######################
MASKED_LTR5="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR5_withprimer.fa"
MASKED_LTR3="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR3_withprimer.fa"

cat "${REF_DIR}/${REF_NAME}" "${REF_DIR}/${VIRUS_NAME}_endU3RU5_withprimer.fa" > "$MASKED_LTR3"
cat "${REF_DIR}/${REF_NAME}" "${REF_DIR}/${VIRUS_NAME}_startU3_withprimer.fa" > "$MASKED_LTR5"

#######################
# Step 2: Index hybrid references
#######################
MMI_LTR3="${MASKED_LTR3%.fa}.mmi"
MMI_LTR5="${MASKED_LTR5%.fa}.mmi"

minimap2 -d "$MMI_LTR3" "$MASKED_LTR3"
minimap2 -d "$MMI_LTR5" "$MASKED_LTR5"

#######################
# Step 3: Map reads for each LTR
#######################
PAF_DIR="${OUT_DIR}/paf"
mkdir -p "$PAF_DIR"

for LTR_NUM in 3 5; do
    if [[ "$LTR_NUM" -eq 3 ]]; then
        MMI_FILE="$MMI_LTR3"
    else
        MMI_FILE="$MMI_LTR5"
    fi

    FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_filtered_size_SUP.fastq"
    SAM_FILE="$(mktemp)"
    BAM_FILE="$(mktemp)"

    PAF_FILE="${PAF_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.paf"

    echo "Mapping $FASTQ_FILE to LTR${LTR_NUM} reference..."
    minimap2 -ax map-ont -t 8 "$MMI_FILE" "$FASTQ_FILE" > "$SAM_FILE"

    # FINAL OUTPUT
    paftools.js sam2paf "$SAM_FILE" > "$PAF_FILE"

    # Cleanup per-LTR
    samtools view -Sb "$SAM_FILE" > "$BAM_FILE"
    rm -f "$SAM_FILE" "$BAM_FILE"
done

#######################
# Cleanup temporary FASTA
#######################
rm -f "$MASKED_LTR3" "$MASKED_LTR5"
