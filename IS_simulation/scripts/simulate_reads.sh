#!/bin/bash
set -euo pipefail

##############################
# Help
##############################
show_help() {
    echo "Usage: $0 -r REF_FASTA -n NUM_SITES -f FLANK_SIZE -l LTR_LINKER_FASTA -d OUT_PREFIX -o OUT_DIR"
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
        o) OUT_DIR="$OPTARG" ;;
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
# Output files (final)
##############################
BED_SORTED="${OUT_DIR}/${OUT_PREFIX}_sorted.bed"
FASTQ_READS="${OUT_DIR}/${OUT_PREFIX}_reads.fq"

##############################
# Temporary files
##############################
BED_RAW="$(mktemp)"
BED_FLANK="$(mktemp)"
FASTA_FLANK="$(mktemp)"
FASTA_READS="$(mktemp)"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

##############################
# Step 1: Generate random integration sites
##############################
python3 "$SCRIPT_DIR/random_sites.py" \
    --fasta "$REF_FASTA" \
    --output "$BED_RAW" \
    --num-sites "$NUM_SITES"

##############################
# Step 2: Sort BED (FINAL BED OUTPUT)
##############################
sort -k1,1 -k2,2n "$BED_RAW" > "$BED_SORTED"

##############################
# Step 3: Generate flanking regions
##############################
python3 "$SCRIPT_DIR/generate_flank_bed.py" \
    "$BED_SORTED" \
    "$BED_FLANK" \
    "$FLANK_SIZE"

##############################
# Step 4: Extract flanking FASTA
##############################
bedtools getfasta \
    -fi "$REF_FASTA" \
    -bed "$BED_FLANK" \
    -fo "$FASTA_FLANK" \
    -nameOnly

##############################
# Step 5: Build simulated reads
##############################
python3 "$SCRIPT_DIR/build_final_random_reads.py" \
    "$FASTA_FLANK" \
    "$LTR_LINKER_FASTA" \
    "$FASTA_READS"

##############################
# Step 6: Convert FASTA to FASTQ (FINAL FASTQ OUTPUT)
##############################
seqtk seq -F 'I' "$FASTA_READS" > "$FASTQ_READS"

##############################
# Cleanup intermediate files
##############################
rm -f \
    "$BED_RAW" \
    "$BED_FLANK" \
    "$FASTA_FLANK" \
    "$FASTA_READS"

echo "DONE"
echo "Final outputs:"
echo "  BED:   $BED_SORTED"
echo "  FASTQ: $FASTQ_READS"
