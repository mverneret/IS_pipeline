#!/bin/bash
set -euo pipefail

##############################
# Help
##############################
show_help() {
    echo "Usage: $0 -r REF_FASTA -n NUM_SITES -f FLANK_SIZE -l LTR_LINKER_FASTA -d OUT_PREFIX -o OUT_DIR"
    echo ""
    echo "Arguments:"
    echo "  -r  Reference genome FASTA (after masking and removing scaffolds)"
    echo "  -n  Number of integration sites to generate"
    echo "  -f  Flanking genome size (bp)"
    echo "  -l  FASTA containing linker + LTR sequences"
    echo "  -d  Output prefix"
    echo "  -o  Output directory"
    echo "  -h  Show this help message"
}

##############################
# Parse arguments
##############################
while getopts "r:n:f:l:d:o:h" opt; do
    case $opt in
        r) REF_FASTA="$OPTARG" ;;
        n) NUM_SITES="$OPTARG" ;;
        f) FLANK_SIZE="$OPTARG" ;;
        l) LTR_LINKER_FASTA="$OPTARG" ;;
        d) OUT_PREFIX="$OPTARG" ;;
        o) OUT_DIR="$OPTARG";;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

##############################
# Check required args
##############################
if [[ -z "${REF_FASTA:-}" || -z "${NUM_SITES:-}" || -z "${FLANK_SIZE:-}" || -z "${LTR_LINKER_FASTA:-}" || -z "${OUT_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$OUT_DIR"

##############################
# Output files
##############################
BED_RAW="${OUT_DIR}/${OUT_PREFIX}.bed"
BED_SORTED="${OUT_DIR}/${OUT_PREFIX}_sorted.bed"
BED_FLANK="${OUT_DIR}/${OUT_PREFIX}_flank.bed"
FASTA_FLANK="${OUT_DIR}/${OUT_PREFIX}_flank.fasta"
FASTA_READS="${OUT_DIR}/${OUT_PREFIX}_reads.fasta"
FASTQ_READS="${OUT_DIR}/${OUT_PREFIX}_reads.fq"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

##############################
# Step 1: Generate random integration sites
##############################
echo "#### Step 1: Generating random integration sites"
python3 $SCRIPT_DIR/random_sites.py \
    --fasta "$REF_FASTA" \
    --output "$BED_RAW" \
    --num-sites "$NUM_SITES"

##############################
# Step 2: Sort BED
##############################
echo "#### Step 2: Sorting BED file"
sort -k1,1 -k2,2n "$BED_RAW" > "$BED_SORTED"

##############################
# Step 3: Count sites per chromosome
##############################
echo "#### Integration sites per chromosome"
cut -f1 "$BED_SORTED" | sort | uniq -c | sort -k2

##############################
# Step 4: Generate flanking regions
##############################
echo "#### Step 4: Generating flanking regions (${FLANK_SIZE} bp)"
python3 $SCRIPT_DIR/generate_flank_bed.py \
    "$BED_SORTED" \
    "$BED_FLANK" \
    "$FLANK_SIZE"

##############################
# Step 5: Extract flanking FASTA
##############################
echo "#### Step 5: Extracting flanking FASTA"
bedtools getfasta \
    -fi "$REF_FASTA" \
    -bed "$BED_FLANK" \
    -fo "$FASTA_FLANK" \
    -nameOnly

##############################
# Step 6: Build simulated reads
##############################
echo "#### Step 6: Building simulated reads"
python3 $SCRIPT_DIR/build_final_random_reads.py \
    "$FASTA_FLANK" \
    "$LTR_LINKER_FASTA" \
    "$FASTA_READS"

##############################
# Step 7: Count read categories
##############################
echo "#### Step 7: Read counts"
for type in \
    "5prime_sens_+" \
    "5prime_antisens_+" \
    "5prime_sens_-" \
    "5prime_antisens_-" \
    "3prime_sens_+" \
    "3prime_antisens_+" \
    "3prime_sens_-" \
    "3prime_antisens_-"
do
    echo -n "Number reads for ${type}: "
    grep -c "$type" "$FASTA_READS" || true
done

##############################
# Step 8: Convert fasta to fastq
##############################
echo "#### Step 8: Convert fa to fq"
seqtk seq -F 'I' "$FASTA_READS" > "$FASTQ_READS"

##############################
# Step 9: Cleanup intermediate files
##############################
echo "####### Cleaning up intermediate files..."
find "$OUT_DIR" -type f ! \( \
  -name "*_sorted.bed" -o \
  -name "*_reads.fq" \
\) -delete

echo "####### DONE"
