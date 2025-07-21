# Integration pipeline
Pipeline inspired from PCIP-seq and Insert-seq that use amplification in long of integration sites so chimeric reads including viral and host sequences. In order to do this libraires prepared amplifying integration sites (IS) with LTR specific primers and primers in Nanopore barcodes. The theorical structure of the reads is shown bellow -> add fixed linker sequences including an identifier called UMI with specific structure  

<img src="image/Image2.png" width="50%">

## Workflow
<img src="image/Image1.png" width="50%">

## Dependancies
- dorado
- bowtie2 (2-2.1.0)
- Nanofilt
- samtools (1.16.1)
- bedtools (v2.30.0)
- python (3.10.14)
- minimap2 (2.26-r1175)
- seqkt (1.5-r133)
- R (> 4.3.1)

## 1- Basecalling
Basecalling was performed on the .fast5 files generated from the different sequencing runs using the Dorado basecaller using super accurate basecalling with GPU acceleration, converting them to .bam formats.
```sh
dorado basecaller --no-trim dna_r10.4.1_e8.2_400bps_sup@v5.0.0 $DATA_DIR > $OUT_BAM_DIR/RUN_${i}_SUP.bam
dorado summary $OUT_BAM_DIR/RUN_${i}_SUP.bam > $OUT_BAM_DIR/RUN_${i}_SUP_summary.tsv
```

## 2- Demultiplexing
Demultiplexing was conducted also using Dorado to separate reads by barcode and remove the Nanopore barcodes at both reads sides.
**#manque étape de merge entre les bam produit par basecalling et le demutiplexing + manque étape bam to fastq after demultiplexing**
```sh
dorado demux --kit-name SQK-NBD114-96 $OUT_BAM_DIR/merged.bam --output-dir $OUTPUT_DIR/demux_both_end --barcode-both-ends --emit-fastq
```

## 3- Filtering steps
Reads were then filtered according to different criteria including :
- Quality > Q20
- Mapping on the LTR5 start or LTR3 end
- Not mapping on provirus without LTR sequences
- Host genome size > 50bp

All these steps can be performed using the ```filter_reads.sh``` script.

**Usage**
```sh
./filter_reads.sh -r <string> -o <string> -f <string> -n <string/int> -q <int> -g <int> -a <int> -l <int> -m <int> -b <int>

Options:
  -r   | Path to LTR_ref_sequences directory
  -o   | Path to output directory
  -f   | Path to demultiplexed FASTQ files directory 
  -n   | Virus reference name (used for fasta ref files and index names 
  -q   | Minimum quality for Nanofilt
  -g   | Minimum genome length in reads
  -a   | Adaptor length
  -l   | LTR5 length in reads (with primer)
  -m   | LTR3 length in reads (with primer)
  -b   | Number of barcodes
```

Outputs:

## 4- Extract UMI 
For the next steps
(inspired from Ivančić et al., 2022)

## 5- Mapping

## 6- Integration sites extraction

___
## Correspondance

___
## Citation

___
