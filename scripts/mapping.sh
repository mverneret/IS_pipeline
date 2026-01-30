#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
    echo "Usage: $0 -r REF_DIR -f REF_NAME -o OUT_DIR -q FASTQ_DIR -n VIRUS_NAME -i SAMPLE_PREFIX [-t TRUE|FALSE]"
}

#######################
# Defaults
#######################
DO_TRIM="FALSE"

#######################
# Parse arguments
#######################
while getopts "r:f:o:q:n:i:t:h" opt; do
    case $opt in
        r) REF_DIR="$OPTARG" ;;
        f) REF_NAME="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        q) FASTQ_DIR="$OPTARG" ;;
        n) VIRUS_NAME="$OPTARG" ;;
        i) SAMPLE_PREFIX="$OPTARG" ;;
        t) DO_TRIM="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

#######################
# Check required arguments
#######################
if [[ -z "${REF_DIR:-}" || -z "${REF_NAME:-}" || -z "${OUT_DIR:-}" || \
      -z "${FASTQ_DIR:-}" || -z "${VIRUS_NAME:-}" || -z "${SAMPLE_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# Step 1: Add LTR to reference
#######################
echo "####### Adding LTR sequences to the host reference"

MASKED_LTR5="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR5_withprimer.fa"
MASKED_LTR3="${REF_DIR}/${REF_NAME%.fa}_${VIRUS_NAME}_LTR3_withprimer.fa"

if [[ ! -f "$MASKED_LTR3" ]]; then
    echo "Creating $MASKED_LTR3"
    cat "${REF_DIR}/${REF_NAME}" \
        "${REF_DIR}/${VIRUS_NAME}_endU3RU5_withprimer.fa" \
        > "$MASKED_LTR3"
else
    echo "Found existing $MASKED_LTR3 — skipping"
fi

if [[ ! -f "$MASKED_LTR5" ]]; then
    echo "Creating $MASKED_LTR5"
    cat "${REF_DIR}/${REF_NAME}" \
        "${REF_DIR}/${VIRUS_NAME}_startU3_withprimer.fa" \
        > "$MASKED_LTR5"
else
    echo "Found existing $MASKED_LTR5 — skipping"
fi

#######################
# Step 2: Index hybrid references
#######################
echo "####### Indexing hybrid references"

MMI_LTR3="${MASKED_LTR3%.fa}.mmi"
MMI_LTR5="${MASKED_LTR5%.fa}.mmi"

if [[ ! -f "$MMI_LTR3" ]]; then
    minimap2 -d "$MMI_LTR3" "$MASKED_LTR3"
else
    echo "Found existing $MMI_LTR3 — skipping"
fi

if [[ ! -f "$MMI_LTR5" ]]; then
    minimap2 -d "$MMI_LTR5" "$MASKED_LTR5"
else
    echo "Found existing $MMI_LTR5 — skipping"
fi

#######################
# Step 2.5: Optional trimming
#######################
if [[ "$DO_TRIM" == "TRUE" ]]; then
    echo "####### Trimming reads with cutadapt"

    for LTR_NUM in 3 5; do
        IN_FASTQ="${FASTQ_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_filtered_size_SUP.fastq"
        OUT_FASTQ="${FASTQ_DIR}/${SAMPLE_PREFIX}_trim_LTR${LTR_NUM}_filtered_size_SUP.fastq"

        if [[ ! -f "$OUT_FASTQ" ]]; then
            echo "Trimming LTR${LTR_NUM}"
            cutadapt \
                -g "GTCTCGTGGGCTCGGTTTVVVVTTVVVVTTVVVVTTVVVVTTTCGCTCTTCCGATCT" \
                -a "AGATCGGAAGAGCGAAABBBBAABBBBAABBBBAABBBBAAACCGAGCCCACGAGAC" \
                -O 15 -e 0.1 \
                -o "$OUT_FASTQ" \
                "$IN_FASTQ"
        else
            echo "Found existing $OUT_FASTQ — skipping"
        fi
    done
else
    echo "####### Trimming disabled"
fi

#######################
# Step 3: Mapping
#######################
for LTR_NUM in 3 5; do
    if [[ "$LTR_NUM" -eq 3 ]]; then
        MMI_FILE="$MMI_LTR3"
    else
        MMI_FILE="$MMI_LTR5"
    fi

    if [[ "$DO_TRIM" == "TRUE" ]]; then
        FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_trim_LTR${LTR_NUM}_filtered_size_SUP.fastq"
    else
        FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_filtered_size_SUP.fastq"
    fi

    SAM_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.sam"
    PAF_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.paf"

    echo "Mapping $FASTQ_FILE to LTR${LTR_NUM}"
    minimap2 -ax map-ont -t 8 "$MMI_FILE" "$FASTQ_FILE" > "$SAM_FILE"

    paftools.js sam2paf "$SAM_FILE" > "$PAF_FILE"
done

#######################
# Step 4: Cleanup
#######################
echo "####### Cleaning up intermediate files"
find "$OUT_DIR" -type f ! -name "*.paf" -delete

echo "####### DONE — $SAMPLE_PREFIX"
