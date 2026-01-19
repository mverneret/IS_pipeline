#!/bin/bash
set -euo pipefail

#######################
# Help function
#######################
show_help() {
    echo "Usage: $0 -s DATA_SAM -o DIR_UMI -a ADAPTER_LENGTH -e MAX_ERROR -i SAMPLE_PREFIX"
}

#######################
# Parse arguments
#######################
while getopts "s:o:a:e:i:h" opt; do
    case $opt in
        s) DATA_SAM="$OPTARG" ;;
        o) DIR_UMI="$OPTARG" ;;
        a) ADAPTER_LENGTH="$OPTARG" ;;
        e) MAX_ERROR="$OPTARG" ;;
        i) SAMPLE_PREFIX="$OPTARG" ;;
        h) show_help; exit 0 ;;
        *) show_help; exit 1 ;;
    esac
done

#######################
# Check required arguments
#######################
if [[ -z "${DATA_SAM:-}" || -z "${DIR_UMI:-}" || -z "${ADAPTER_LENGTH:-}" || -z "${MAX_ERROR:-}" || -z "${SAMPLE_PREFIX:-}" ]]; then
    echo "Error: Missing required arguments"
    show_help
    exit 1
fi

mkdir -p "$DIR_UMI"

#######################
# Script directory
#######################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#######################
# 1. Separate reads by LTR and orientation
#######################
for a in "endU3RU5" "startU3"; do
    if [[ $a == "endU3RU5" ]]; then
        LTR="LTR3"
    else
        LTR="LTR5"
    fi

    FORWARD_OUT="${DIR_UMI}/${SAMPLE_PREFIX}_${LTR}_SUP_fwd.fasta"
    REVERSE_OUT="${DIR_UMI}/${SAMPLE_PREFIX}_${LTR}_SUP_rev.fasta"

    echo "Processing $LTR reads (forward and reverse)..."

    samtools view -F 20 "${DATA_SAM}/${SAMPLE_PREFIX}_mapping_${a}_SUP.sam" \
        | grep -v '^@' \
        | awk '{print ">" $1 "\n" $10}' > "$FORWARD_OUT"

    samtools view -f 16 -F 4 "${DATA_SAM}/${SAMPLE_PREFIX}_mapping_${a}_SUP.sam" \
        | grep -v '^@' \
        | awk '{print ">" $1 "\n" $10}' > "$REVERSE_OUT"
done

#######################
# 2. Run UMI extraction
#######################
for LTR_NUM in 3 5; do
    FWD_FILE="${DIR_UMI}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_SUP_fwd.fasta"
    REV_FILE="${DIR_UMI}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_SUP_rev.fasta"
    OUT_FILE="${DIR_UMI}/${SAMPLE_PREFIX}_LTR${LTR_NUM}_UMI.fasta"

    echo "Extracting UMI for LTR${LTR_NUM}..."
    python3 "${SCRIPT_DIR}/insert_seq_extract_umi_modif.py" \
        --max-error "$MAX_ERROR" \
        --adapter-length "$ADAPTER_LENGTH" \
        -o "$OUT_FILE" \
        -ifwd "$FWD_FILE" \
        -irev "$REV_FILE"
done

#######################
# Cleanup
#######################
rm -f \
    "${DIR_UMI}/${SAMPLE_PREFIX}_LTR3_SUP_fwd.fasta" \
    "${DIR_UMI}/${SAMPLE_PREFIX}_LTR3_SUP_rev.fasta" \
    "${DIR_UMI}/${SAMPLE_PREFIX}_LTR5_SUP_fwd.fasta" \
    "${DIR_UMI}/${SAMPLE_PREFIX}_LTR5_SUP_rev.fasta"

echo "DONE"
echo "Final outputs:"
echo "  ${DIR_UMI}/${SAMPLE_PREFIX}_LTR3_UMI.fasta"
echo "  ${DIR_UMI}/${SAMPLE_PREFIX}_LTR5_UMI.fasta"
