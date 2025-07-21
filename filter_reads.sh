#!/bin/bash

# Exit on error
set -euo pipefail

# Help function
show_help() {
  echo "Usage: $0 -r REF_DIR -o OUT_DIR -f FASTQ_DEMUX -n REF_NAME -q MIN_QUALITY -g LENGTH_GENOME -a LENGTH_ADAPTOR -l LENGTH_LTR5 -m LENGTH_LTR3 -b NB_BARCODES"
  echo ""
  echo "Arguments:"
  echo "  -r  Path to reference directory"
  echo "  -o  Output directory"
  echo "  -f  Path to demultiplexed FASTQ files"
  echo "  -n  Reference name (used for fasta and index names)"
  echo "  -q  Minimum quality for NanoFilt"
  echo "  -g  Min genome length in the reads"
  echo "  -a  Adaptor length"
  echo "  -l  LTR5 length in the reads (with primers)"
  echo "  -m  LTR3 length in the reads (with primers)"
  echo "  -b  Number of barcodes (e.g., 12 for barcode01 to barcode12)"
  echo "  -h  Show this help message"
}

# Parse options
while getopts "r:o:f:n:q:g:a:l:m:b:h" opt; do
  case $opt in
    r) REF_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    f) FASTQ_DEMUX="$OPTARG" ;;
    n) ref_name="$OPTARG" ;;
    q) min_quality="$OPTARG" ;;
    g) length_genome="$OPTARG" ;;
    a) length_adaptor="$OPTARG" ;;
    l) length_LTR5="$OPTARG" ;;
    m) length_LTR3="$OPTARG" ;;
    b) nb_barcodes="$OPTARG" ;;
    h) show_help; exit 0 ;;
    *) show_help; exit 1 ;;
  esac
done

# Check all required args
if [[ -z "${REF_DIR:-}" || -z "${OUT_DIR:-}" || -z "${FASTQ_DEMUX:-}" || -z "${ref_name:-}" || -z "${min_quality:-}" || -z "${length_genome:-}" || -z "${length_adaptor:-}" || -z "${length_LTR5:-}" || -z "${length_LTR3:-}" || -z "${nb_barcodes:-}" ]]; then
  echo "Error: Missing required arguments"
  show_help
  exit 1
fi

# Build index files if they don't exist
echo "####### 0. Index LTR and provirus references"

for ref in "endU3RU5" "startU3" "provirus_wo_LTR"; do
  index="${REF_DIR}/${ref_name}_${ref}_index.1.bt2"
  fasta="${REF_DIR}/${ref_name}_${ref}.fasta"
  if [ ! -f "$index" ]; then
    echo "Index for ${fasta} does not exist. Building index..."
    bowtie2-build "$fasta" "${fasta%.fasta}_index"
  else
    echo "Index for ${fasta} already exists. Skipping."
  fi
done

# Process each barcode
for i in $(seq -w 1 $nb_barcodes); do
  echo "####### 1. Filter read quality - barcode$i"
  NanoFilt -q "$min_quality" "${FASTQ_DEMUX}/SQK-NBD114-96_barcode${i}.fastq" > "${FASTQ_DEMUX}/SQK-NBD114-96_barcode${i}_Q${min_quality}.fastq"

  echo "####### 2. Mapping reads on LTRs - barcode$i"
  # LTR3
  bowtie2 --sensitive-local -x "${REF_DIR}/${ref_name}_endU3RU5_index" -U "${FASTQ_DEMUX}/SQK-NBD114-96_barcode${i}_Q${min_quality}.fastq" -S "${OUT_DIR}/barcode${i}_mapping_endU3RU5_SUP.sam" -N 1
  samtools view -h -F 4 "${OUT_DIR}/barcode${i}_mapping_endU3RU5_SUP.sam" > "${OUT_DIR}/barcode${i}_mapped_endU3RU5_SUP.sam"
  samtools fastq "${OUT_DIR}/barcode${i}_mapped_endU3RU5_SUP.sam" > "${OUT_DIR}/barcode${i}_mapped_endU3RU5_SUP.fastq"
  # LTR5
  bowtie2 --sensitive-local -x "${REF_DIR}/${ref_name}_startU3_index" -U "${FASTQ_DEMUX}/SQK-NBD114-96_barcode${i}_Q${min_quality}.fastq" -S "${OUT_DIR}/barcode${i}_mapping_startU3_SUP.sam" -N 1
  samtools view -h -F 4 "${OUT_DIR}/barcode${i}_mapping_startU3_SUP.sam" > "${OUT_DIR}/barcode${i}_mapped_startU3_SUP.sam"
  samtools fastq "${OUT_DIR}/barcode${i}_mapped_startU3_SUP.sam" > "${OUT_DIR}/barcode${i}_mapped_startU3_SUP.fastq"

  echo "####### 3. Mapping reads on provirus w/o LTR - barcode$i"
  # LTR3
  bowtie2 --sensitive-local -x "${REF_DIR}/${ref_name}_provirus_wo_LTR_index" -U "${OUT_DIR}/barcode${i}_mapped_endU3RU5_SUP.fastq" -S "${OUT_DIR}/barcode${i}_mapped_endU3RU5_mapping_${ref_name}_provirus_wo_LTR_SUP.sam" -N 1
  samtools view -h -f 4 "${OUT_DIR}/barcode${i}_mapped_endU3RU5_mapping_${ref_name}_provirus_wo_LTR_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR3_SUP.sam"
  # LTR5
  bowtie2 --sensitive-local -x "${REF_DIR}/${ref_name}_provirus_wo_LTR_index" -U "${OUT_DIR}/barcode${i}_mapped_startU3_SUP.fastq" -S "${OUT_DIR}/barcode${i}_mapped_startU3_mapping_${ref_name}_provirus_wo_LTR_SUP.sam" -N 1
  samtools view -h -f 4 "${OUT_DIR}/barcode${i}_mapped_startU3_mapping_${ref_name}_provirus_wo_LTR_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR5_SUP.sam"

  echo "####### 4. Filter read size - barcode$i"
  # LTR3
  min_length_LTR3=$((length_genome + length_adaptor + length_LTR3))
  awk -v min_len="$min_length_LTR3" 'length($10) >= min_len || $1 ~ /^@/' "${OUT_DIR}/barcode${i}_LTR3_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR3_filtered_size_SUP.sam"
  samtools fastq "${OUT_DIR}/barcode${i}_LTR3_filtered_size_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR3_filtered_size_SUP.fastq"
  # LTR5
  min_length_LTR5=$((length_genome + length_adaptor + length_LTR5))
  awk -v min_len="$min_length_LTR5" 'length($10) >= min_len || $1 ~ /^@/' "${OUT_DIR}/barcode${i}_LTR5_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR5_filtered_size_SUP.sam"
  samtools fastq "${OUT_DIR}/barcode${i}_LTR5_filtered_size_SUP.sam" > "${OUT_DIR}/barcode${i}_LTR5_filtered_size_SUP.fastq"
done

echo "####### DONE."
