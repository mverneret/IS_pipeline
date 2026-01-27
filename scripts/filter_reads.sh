#!/bin/bash

# Exit on error
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

# Check required args
if [[ -z "${REF_DIR:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_FILE:-}" || -z "${REF_NAME:-}" || -z "${SAMPLE_NAME:-}" || -z "${MIN_QUALITY:-}" || -z "${LENGTH_GENOME:-}" || -z "${LENGTH_ADAPTOR:-}" || -z "${LENGTH_LTR5:-}" || -z "${LENGTH_LTR3:-}" ]]; then
  echo "Error: Missing required arguments"
  show_help
  exit 1
fi

mkdir -p "$OUT_DIR"

#######################
# 0. Build reference indexes if missing
#######################
echo "####### 0. Index LTR and provirus references"
for ref in "endU3RU5" "startU3" "provirus_wo_LTR"; do
  INDEX="${REF_DIR}/${REF_NAME}_${ref}_index.1.bt2"
  FASTA="${REF_DIR}/${REF_NAME}_${ref}.fasta"
  if [ ! -f "$INDEX" ]; then
    echo "Index for ${FASTA} does not exist. Building index..."
    bowtie2-build "$FASTA" "${FASTA%.fasta}_index"
  else
    echo "Index for ${FASTA} already exists. Skipping."
  fi
done

#######################
# 1. Filter read quality
#######################
echo "####### 1. Filter read quality - ${SAMPLE_NAME}"
FILTERED_FASTQ="${OUT_DIR}/${SAMPLE_NAME}_Q${MIN_QUALITY}.fastq"
NanoFilt -q "$MIN_QUALITY" "$FASTQ_FILE" > "$FILTERED_FASTQ"

#######################
# 2. Mapping reads on LTRs
#######################
echo "####### 2. Mapping reads on LTRs - ${SAMPLE_NAME}"

# LTR3 (endU3RU5)
echo "####### Mapping on LTR3"
SAM_END="${OUT_DIR}/${SAMPLE_NAME}_mapping_endU3RU5_SUP.sam"
bowtie2 --sensitive-local -x "${REF_DIR}/${REF_NAME}_endU3RU5_index" -U "$FILTERED_FASTQ" -S "$SAM_END" -N 1
samtools view -h -F 4 "$SAM_END" > "${OUT_DIR}/${SAMPLE_NAME}_mapped_endU3RU5_SUP.sam"
samtools fastq "${OUT_DIR}/${SAMPLE_NAME}_mapped_endU3RU5_SUP.sam" > "${OUT_DIR}/${SAMPLE_NAME}_mapped_endU3RU5_SUP.fastq"

# LTR5 (startU3)
echo "####### Mapping on LTR5"
SAM_START="${OUT_DIR}/${SAMPLE_NAME}_mapping_startU3_SUP.sam"
bowtie2 --sensitive-local -x "${REF_DIR}/${REF_NAME}_startU3_index" -U "$FILTERED_FASTQ" -S "$SAM_START" -N 1
samtools view -h -F 4 "$SAM_START" > "${OUT_DIR}/${SAMPLE_NAME}_mapped_startU3_SUP.sam"
samtools fastq "${OUT_DIR}/${SAMPLE_NAME}_mapped_startU3_SUP.sam" > "${OUT_DIR}/${SAMPLE_NAME}_mapped_startU3_SUP.fastq"

#######################
# 3. Mapping reads on provirus without LTR
#######################
echo "####### 3. Mapping reads on provirus without LTR - ${SAMPLE_NAME}"

SAM_PRO_END="${OUT_DIR}/${SAMPLE_NAME}_mapped_startU3_mapping_${REF_NAME}_provirus_wo_LTR_SUP.sam"
bowtie2 --sensitive-local -x "${REF_DIR}/${REF_NAME}_provirus_wo_LTR_index" \
  -U "${OUT_DIR}/${SAMPLE_NAME}_mapped_endU3RU5_SUP.fastq" -S "$SAM_PRO_END" -N 1
samtools view -h -f 4 "$SAM_PRO_END" > "${OUT_DIR}/${SAMPLE_NAME}_LTR3_SUP.sam"

SAM_PRO_START="${OUT_DIR}/${SAMPLE_NAME}_mapped_endU3RU5_mapping_${REF_NAME}_provirus_wo_LTR_SUP.sam"
bowtie2 --sensitive-local -x "${REF_DIR}/${REF_NAME}_provirus_wo_LTR_index" \
  -U "${OUT_DIR}/${SAMPLE_NAME}_mapped_startU3_SUP.fastq" -S "$SAM_PRO_START" -N 1
samtools view -h -f 4 "$SAM_PRO_START" > "${OUT_DIR}/${SAMPLE_NAME}_LTR5_SUP.sam"

#######################
# 4. Filter read size
#######################
echo "####### 4. Filter read size - ${SAMPLE_NAME}"

MIN_LEN_LTR3=$((LENGTH_GENOME + LENGTH_ADAPTOR + LENGTH_LTR3))
awk -v min_len="$MIN_LEN_LTR3" 'length($10) >= min_len || $1 ~ /^@/' "${OUT_DIR}/${SAMPLE_NAME}_LTR3_SUP.sam" \
  > "${OUT_DIR}/${SAMPLE_NAME}_LTR3_filtered_size_SUP.sam"
samtools fastq "${OUT_DIR}/${SAMPLE_NAME}_LTR3_filtered_size_SUP.sam" > "${OUT_DIR}/${SAMPLE_NAME}_LTR3_filtered_size_SUP.fastq"

MIN_LEN_LTR5=$((LENGTH_GENOME + LENGTH_ADAPTOR + LENGTH_LTR5))
awk -v min_len="$MIN_LEN_LTR5" 'length($10) >= min_len || $1 ~ /^@/' "${OUT_DIR}/${SAMPLE_NAME}_LTR5_SUP.sam" \
  > "${OUT_DIR}/${SAMPLE_NAME}_LTR5_filtered_size_SUP.sam"
samtools fastq "${OUT_DIR}/${SAMPLE_NAME}_LTR5_filtered_size_SUP.sam" > "${OUT_DIR}/${SAMPLE_NAME}_LTR5_filtered_size_SUP.fastq"

#######################
# 5. Cleanup intermediate files
#######################
echo "####### Cleaning up intermediate files..."

# Keep only mapping SAMs and filtered FASTQ
find "$OUT_DIR" -type f ! \( \
  -name "*_mapping_endU3RU5_SUP.sam" -o \
  -name "*_mapping_startU3_SUP.sam" -o \
  -name "*_LTR3_filtered_size_SUP.fastq" -o \
  -name "*_LTR5_filtered_size_SUP.fastq" \
\) -delete

echo "####### DONE - ${SAMPLE_NAME}"
