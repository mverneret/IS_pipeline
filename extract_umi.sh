#!/bin/bash
set -euo pipefail

# Check and parse arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <DATA_SAM> <DIR_UMI> <nb_barcodes> <adapter_length> <max_error>"
    exit 1
fi

DATA_SAM="$1"
DIR_UMI="$2"
nb_barcodes="$3"
adapter_length="$4"
max_error="$5"

# Step 1: Separate reads by LTR and orientation
for i in $(seq -w 1 $nb_barcodes); do
    for a in endU3RU5 startU3; do
        if [[ $a == "endU3RU5" ]]; then
            LTR="LTR3"
        elif [[ $a == "startU3" ]]; then
            LTR="LTR5"
        fi
        samtools view -F 20 ${DATA_SAM}/barcode${i}_mapping_${a}_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > ${DIR_UMI}/barcode${i}_${LTR}_SUP_fwd.fasta
        samtools view -f 16 -F 4 ${DATA_SAM}/barcode${i}_mapping_${a}_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > ${DIR_UMI}/barcode${i}_${LTR}_SUP_rev.fasta
    done
done

# Step 2: Run UMI extraction for each barcode and LTR
for i in $(seq -w 1 $nb_barcodes); do
    for a in {3,5}; do
        python3 insert_seq_extract_umi_modif.py \
            --max-error "$max_error" \
            --adapter-length "$adapter_length" \
            -o "${DIR_UMI}/barcode${i}_LTR${a}_UMI.fasta" \
            -ifwd "${DIR_UMI}/barcode${i}_LTR${a}_SUP_fwd.fasta" \
            -irev "${DIR_UMI}/barcode${i}_LTR${a}_SUP_rev.fasta"
    done
done
