#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
    echo "Usage: $0 -r REF_DIR -f REF_NAME -t TE_ANNOT -o OUT_DIR -q FASTQ_DIR -c PREFIX_CHR -n VIRUS_NAME -i SAMPLE_PREFIX"
    echo ""
    echo "Arguments:"
    echo "  -r  Reference directory"
    echo "  -f  Reference host genome FASTA file name"
    echo "  -t  TE annotation BED file"
    echo "  -o  Output directory"
    echo "  -q  Directory containing FASTQ files"
    echo "  -c  Prefix of chromosomes to keep"
    echo "  -n  Virus reference name (for LTR files)"
    echo "  -i  Sample prefix for output files"
    echo "  -h  Show this help message"
}

#######################
# Parse arguments
#######################
while getopts "r:f:t:o:q:c:n:i:h" opt; do
    case $opt in
        r) REF_DIR="$OPTARG" ;;
        f) REF_NAME="$OPTARG" ;;
        t) TE_ANNOT="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        q) FASTQ_DIR="$OPTARG" ;;
        c) PREFIX_CHR="$OPTARG" ;;
        n) VIRUS_NAME="$OPTARG" ;;
        i) SAMPLE_PREFIX="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

# Check required arguments
if [[ -z "${REF_DIR:-}" || -z "${REF_NAME:-}" || -z "${TE_ANNOT:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_DIR:-}" || -z "${PREFIX_CHR:-}" || -z "${VIRUS_NAME:-}" || -z "${SAMPLE_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# Step 1: Mask reference genome
#######################
echo "####### Mask host reference genome with TE annotation"
MASKED_REF="${REF_DIR}/${REF_NAME%.fa}_masked.fa"
bedtools maskfasta -fi "${REF_DIR}/${REF_NAME}" -bed "${REF_DIR}/${TE_ANNOT}" -fo "$MASKED_REF"

#######################
# Step 2: Remove scaffolds
#######################
echo "####### Remove scaffolds from the host reference genome"
NOSCAFF_REF="${REF_DIR}/${REF_NAME%.fa}_noscaffold_masked.fa"
awk -v prefix_chr="^>${PREFIX_CHR}" '
BEGIN { keep=0 }
/^>/ { keep = ($0 ~ prefix_chr) }
keep { print }
' "$MASKED_REF" > "$NOSCAFF_REF"

#######################
# Step 3: Add LTR to reference
#######################
echo "####### Add LTR sequences to the host reference genome"
MASKED_U5="${REF_DIR}/${REF_NAME%.fa}_masked_U5_withprimer.fa"
MASKED_U3="${REF_DIR}/${REF_NAME%.fa}_masked_U3_withprimer.fa"
cat "$NOSCAFF_REF" "${REF_DIR}/${VIRUS_NAME}_endU3RU5_withprimer.fa" > "$MASKED_U5"
cat "$NOSCAFF_REF" "${REF_DIR}/${VIRUS_NAME}_startU3_withprimer.fa" > "$MASKED_U3"

#######################
# Step 4: Index hybrid references
#######################
minimap2 -d "${MASKED_U3%.fa}.mmi" "$MASKED_U3"
minimap2 -d "${MASKED_U5%.fa}.mmi" "$MASKED_U5"

#######################
# Step 5: Map reads for each LTR
#######################
for LTR_NUM in 3 5; do
    if [ "$LTR_NUM" -eq 3 ]; then
        MMI_FILE="${MASKED_U5%.fa}.mmi"
    else
        MMI_FILE="${MASKED_U3%.fa}.mmi"
    fi

    FASTQ_FILE="${FASTQ_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_filtered_size_SUP.fastq"
    SAM_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.sam"
    BAM_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.bam"
    SORTED_BAM_FILE="${OUT_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_sorted_SUP.bam"

    echo "####### Mapping $FASTQ_FILE to LTR${LTR_NUM} reference..."
    minimap2 -ax map-ont -t 8 "$MMI_FILE" "$FASTQ_FILE" > "$SAM_FILE"

    # Convert to PAF
    PAF_DIR="${OUT_DIR}/paf"
    mkdir -p "$PAF_DIR"
    paftools.js sam2paf "$SAM_FILE" > "${PAF_DIR}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_mapped_${REF_NAME}_SUP.paf"

    # Convert SAM -> BAM -> sorted BAM
    samtools view -@ 4 -Sb "$SAM_FILE" > "$BAM_FILE"
    samtools sort -@ 4 "$BAM_FILE" -o "$SORTED_BAM_FILE"
    samtools index -@ 4 "$SORTED_BAM_FILE"
done

echo "####### DONE - $SAMPLE_PREFIX"
