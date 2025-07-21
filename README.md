# Integration pipeline

Structure des reads thérorique + image pipeline global

## 1- Basecalling
Run on GPU
Model : dna_r10.4.1_e8.2_400bps_sup@v5.0.0
dorado basecaller --no-trim dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $DATA_DIR > $OUT_BAM_DIR/RUN_02_SUP.bam
dorado summary $OUT_BAM_DIR/RUN_02_SUP.bam > $OUT_BAM_DIR/RUN_02_SUP_summary.tsv

## 2- Demultiplexing
#manque étape de merge entre les bam produit par basecalling et le demutiplexing
dorado demux --kit-name SQK-NBD114-96 $OUT_BAM_DIR/merged.bam --output-dir $OUTPUT_DIR/demux_both_end --barcode-both-ends --emit-fastq

## 3- Filtering steps
./filter_reads.sh \
  -r /path/to/ref \
  -o /path/to/output \
  -f /path/to/fastq_demux \
  -n 2369 \
  -q 20 \
  -g 50 \
  -a 57 \
  -l 37 \
  -m 162 \
  -b 20

## 4- Extract UMI 
(inspired from Ivančić et al., 2022)

## 5- Mapping

## 6- Integration sites extraction

