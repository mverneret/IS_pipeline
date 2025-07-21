#!/bin/bash

show_help() {
  echo "Usage: $0 <REF_DIR> <REF_name> <TE_annot> <OUT_DIR> <FASTQ_DIR> <prefix_chr> <virus_name> <nb_barcodes>"
}

REF_DIR=$1
REF=$2
TE_annot=$3
OUT_DIR=$4
FASTQ=$5
prefix_chr=$6
virus_name=$7
nb_barcodes=$8

# Create output directory if needed
mkdir -p "$OUT_DIR"

# Step 1: Mask reference genome
bedtools maskfasta -fi "${REF_DIR}/${REF}" -bed "${REF_DIR}/${TE_annot}" -fo "${REF_DIR}/${REF%.fa}_masked.fa"

# Step 2: Remove scaffolds
awk -v prefix_chr="^>${prefix_chr}" '
BEGIN { keep=0 }
/^>/ { keep = ($0 ~ prefix_chr) }
keep { print }
' "${REF_DIR}/${REF%.fa}_masked.fa" > "${REF_DIR}/${REF%.fa}_noscaffold_masked.fa"

# Step 3: Add the LTR to the ref
cat "${REF_DIR}/${REF%.fa}_noscaffold_masked.fa" "${REF_DIR}/${virus_name}_endU3RU5_withprimer.fa" > "${REF_DIR}/${REF%.fa}_masked_U5_withprimer.fa"
cat "${REF_DIR}/${REF%.fa}_noscaffold_masked.fa" "${REF_DIR}/${virus_name}_startU3_withprimer.fa" > "${REF_DIR}/${REF%.fa}_masked_U3_withprimer.fa"

# Step 4: Index hybrid references 
minimap2 -d "${REF_DIR}/${REF%.fa}_masked_U3_withprimer.mmi" "${REF_DIR}/${REF%.fa}_masked_U3_withprimer.fa"
minimap2 -d "${REF_DIR}/${REF%.fa}_masked_U5_withprimer.mmi" "${REF_DIR}/${REF%.fa}_masked_U5_withprimer.fa"

# Step 5: Mapping reads for each barcode
for i in $(seq -w 1 $nb_barcodes); do
  for a in 3 5; do
    if [ "$a" -eq 3 ]; then
      MMI_FILE="${REF_DIR}/${REF%.fa}_masked_U5_withprimer.mmi"
    else
      MMI_FILE="${REF_DIR}/${REF%.fa}_masked_U3_withprimer.mmi"
    fi

    FASTQ_FILE="${FASTQ}/barcode${i}_LTR${a}_filtered_size_SUP.fastq"
    SAM_FILE="${OUT_DIR}/barcode${i}_LTR${a}_mapped_${REF}_SUP.sam"
    BAM_FILE="${OUT_DIR}/barcode${i}_LTR${a}_mapped_${REF}_SUP.bam"
    SORTED_BAM_FILE="${OUT_DIR}/barcode${i}_LTR${a}_mapped_${REF}_sorted_SUP.bam"

    # Align with minimap2
    minimap2 -ax map-ont -t 8 "$MMI_FILE" "$FASTQ_FILE" > "$SAM_FILE"

    # Convert SAM to PAF if needed
    OUT_DIR_PAF="${OUT_DIR}/paf"
    mkdir -p "$OUT_DIR_PAF"
    paftools.js sam2paf "$SAM_FILE" > "${OUT_DIR_PAF}/barcode${i}_LTR${a}_mapped_${REF}_SUP.paf"

    # Convert SAM to BAM and sort/index
    samtools view -@ 4 -Sb "$SAM_FILE" > "$BAM_FILE"
    samtools sort -@ 4 "$BAM_FILE" -o "$SORTED_BAM_FILE"
    samtools index -@ 4 "$SORTED_BAM_FILE"
  done
done
