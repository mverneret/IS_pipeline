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

## 3'- Extract UMI (inspired from Ivančić et al., 2022)
DATA_SAM=/beegfs/project/rer/Seq/2025ONTPROM01/results/bowtie2/JSRV
DIR_UMI=/beegfs/project/rer/Seq/2025ONTPROM01/results/extract_UMI/JSRV
for i in {01..26}; do
    for a in {endU3RU5,startU3}; do

  #Determine the LTR label based on the value of 'a'
        if [[ $a == "endU3RU5" ]]; then
            LTR="LTR3"
        elif [[ $a == "startU3" ]]; then
            LTR="LTR5"
        fi

                samtools view -F 20 ${DATA_SAM}/barcode${i}_mapping_${a}_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > ${DIR_UMI}/barcode${i}_${LTR}_SUP_fwd.fasta
                samtools view -f 16 -F 4 ${DATA_SAM}/barcode${i}_mapping_${a}_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > ${DIR_UMI}/barcode${i}_${LTR}_SUP_rev.fasta
        done
done

#Define parameters
max_error=0
adapter_length=57

#Loop through barcode numbers and LTR types
for i in {01..26}; do
    for a in {3,5}; do

        # Run the Python script with the specified parameters
        python3 insert_seq_extract_umi_modif3.py \
            --max-error $max_error \
            --adapter-length $adapter_length \
            -o ${DIR_UMI}/barcode${i}_LTR${a}_UMI.fasta \
        -ifwd ${DIR_UMI}/barcode${i}_LTR${a}_SUP_fwd.fasta \
        -irev ${DIR_UMI}/barcode${i}_LTR${a}_SUP_rev.fasta

        # Check if the command was successful
        if [[ $? -eq 0 ]]; then
            echo "Successfully processed $input_file and saved to $output_file"
        else
            echo "Error processing $input_file"
        fi
    done
done

