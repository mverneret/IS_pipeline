#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
  echo "Usage: $0 -r REF_DIR -o OUT_DIR -f FASTQ_FILE -n REF_NAME -i SAMPLE_NAME -q MIN_QUALITY -g LENGTH_GENOME -a LENGTH_ADAPTOR -l LENGTH_LTR5 -m LENGTH_LTR3"
}

#######################
# Parse options
#######################
while getopts "r:o:f:n:i:q:g:a:l:m:h" opt; do
  case $opt in
    r) REF_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    f) FASTQ_FILE="$OPTARG" ;;
    n) REF_NAME="$OPTARG" ;;
    i) SAMPLE_NAME="$OPTARG" ;;
    q) MIN_QUALITY="$OPTARG" ;;
    g) LENGTH_GENOME="$OPTARG" ;;
    a) LENGTH_ADAPTOR="$OPTARG" ;;
    l) LENGTH_LTR5="$OPTARG" ;;
    m) LENGTH_LTR3="$OPTARG" ;;
    h) show_help; exit 0 ;;
    *) show_help; exit 1 ;;
  esac
done

#######################
# Check required args
#######################
if [[ -z "${REF_DIR:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_FILE:-}" || -z "${REF_NAME:-}" || -z "${SAMPLE_NAME:-}" || -z "${MIN_QUALITY:-}" || -z "${LENGTH_GENOME:-}" || -z "${LENGTH_ADAPTOR:-}" || -z "${LENGTH_LTR5:-}" || -z "${LENGTH_LTR3:-}" ]]; then
  echo "Error: Missing required arguments"
  show_help
  exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# 0. Build reference indexes if missing
#######################
for ref in "endU3RU5" "startU3" "provirus_wo_LTR"; do
  FASTA="${REF_DIR}/${REF_NAME}_${ref}.fasta"
  INDEX="${FASTA%.fasta}_index.1.bt2"
  if [[ ! -f "$INDEX" ]]; then
    bowtie2-build "$FASTA" "${FASTA%.fasta}_index"
  fi
done

#######################
# 1. Filter read quality 
#######################
FILTERED_FASTQ="$(mktemp)"
NanoFilt -q "$MIN_QUALITY" "$FASTQ_FILE" > "$FILTERED_FASTQ"

#######################
# 2. Mapping reads on LTRs (FINAL SAMs)
#######################

# LTR3
SAM_LTR3="${OUT_DIR}/${SAMPLE_NAME}_mapping_LTR3_SUP.sam"
bowtie2 --sensitive-local \
  -x "${REF_DIR}/${REF_NAME}_endU3RU5_index" \
  -U "$FILTERED_FASTQ" \
  -S "$SAM_LTR3" -N 1

samtools view -h -F 4 "$SAM_LTR3" > "$(mktemp)" && mv "$(mktemp)" "$SAM_LTR3"

# LTR5
SAM_LTR5="${OUT_DIR}/${SAMPLE_NAME}_mapping_LTR5_SUP.sam"
bowtie2 --sensitive-local \
  -x "${REF_DIR}/${REF_NAME}_startU3_index" \
  -U "$FILTERED_FASTQ" \
  -S "$SAM_LTR5" -N 1

samtools view -h -F 4 "$SAM_LTR5" > "$(mktemp)" && mv "$(mktemp)" "$SAM_LTR5"

#######################
# 3. Mapping on provirus
#######################
TMP_LTR3_FASTQ="$(mktemp)"
TMP_LTR5_FASTQ="$(mktemp)"

samtools fastq "$SAM_LTR3" > "$TMP_LTR3_FASTQ"
samtools fastq "$SAM_LTR5" > "$TMP_LTR5_FASTQ"

TMP_LTR3_SAM="$(mktemp)"
TMP_LTR5_SAM="$(mktemp)"

bowtie2 --sensitive-local \
  -x "${REF_DIR}/${REF_NAME}_provirus_wo_LTR_index" \
  -U "$TMP_LTR3_FASTQ" \
  -S "$TMP_LTR3_SAM" -N 1

bowtie2 --sensitive-local \
  -x "${REF_DIR}/${REF_NAME}_provirus_wo_LTR_index" \
  -U "$TMP_LTR5_FASTQ" \
  -S "$TMP_LTR5_SAM" -N 1

#######################
# 4. Filter read size (FINAL FASTQs)
#######################

# LTR3
MIN_LEN_LTR3=$((LENGTH_GENOME + LENGTH_ADAPTOR + LENGTH_LTR3))
OUT_LTR3_FASTQ="${OUT_DIR}/${SAMPLE_NAME}_LTR3_filtered_size_SUP.fastq"

awk -v min_len="$MIN_LEN_LTR3" 'length($10) >= min_len || $1 ~ /^@/' "$TMP_LTR3_SAM" \
  | samtools fastq - > "$OUT_LTR3_FASTQ"

# LTR5
MIN_LEN_LTR5=$((LENGTH_GENOME + LENGTH_ADAPTOR + LENGTH_LTR5))
OUT_LTR5_FASTQ="${OUT_DIR}/${SAMPLE_NAME}_LTR5_filtered_size_SUP.fastq"

awk -v min_len="$MIN_LEN_LTR5" 'length($10) >= min_len || $1 ~ /^@/' "$TMP_LTR5_SAM" \
  | samtools fastq - > "$OUT_LTR5_FASTQ"

#######################
# Cleanup
#######################
rm -f \
  "$FILTERED_FASTQ" \
  "$TMP_LTR3_FASTQ" \
  "$TMP_LTR5_FASTQ" \
  "$TMP_LTR3_SAM" \
  "$TMP_LTR5_SAM"

echo "DONE"
echo "Final outputs:"
echo "  $SAM_LTR3"
echo "  $SAM_LTR5"
echo "  $OUT_LTR3_FASTQ"
echo "  $OUT_LTR5_FASTQ"

