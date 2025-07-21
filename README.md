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
Basecalling was performed on the files generated from the different sequencing runs using the Dorado basecaller using super accurate basecalling with GPU acceleration, converting .fast5 to .bam formats.
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

All these steps can be performed using the ```filter_reads.sh``` script that will loop on the different barcodes and deal with LTR5 and LTR3 reads.

**Usage**
```sh
./filter_reads.sh -r <string> -o <string> -f <string> -n <string/int> -q <int> -g <int> -a <int> -l <int> -m <int> -b <int>

Options:
  -r <REF_DIR>         | Path to LTR_ref_sequences directory
  -o <OUT_DIR>         | Path to output directory
  -f <FASTQ_DEMUX>     | Path to demultiplexed FASTQ files directory 
  -n <ref_name>        | Virus reference name (used for fasta ref files and index names 
  -q <min_quality>     | Minimum quality for Nanofilt
  -g <length_genome>   | Minimum genome length in reads
  -a <length_adaptor>  | Adaptor length
  -l <length_LTR5>     | LTR5 length in reads (with primer)
  -m <length LTR3>     | LTR3 length in reads (with primer)
  -b <nb_barcode>      | Number of barcodes
```

Outputs in ```OUT_DIR``` (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```barcode${i}_mapping_${a}_SUP.sam``` : Mapping results on start LTR5 or end LTR3
- ```barcode${i}_mapped_${a}_SUP.sam|fastq``` : Only mapped reads on start LTR5 or end LTR3
- ```barcode${i}_mapped_${a}_mapping_${ref_name}_provirus_wo_LTR_SUP.sam``` : Mapping results on provirus w/o LTR sequences
- ```barcode${i}_${a}_SUP.sam``` : Only non-mapped reads on provirus w/o LTR sequences
- ```barcode${i}_${a}_filtered_size_SUP.sam|fastq``` : Reads after all the steps + size filtering

## 4- Extract UMI 
In order to remove PCR duplicates for clonality quantification in the Step 6- it was necessary to extract UMI sequences from all the reads. We thus modified a python script from the INSERT-seq pipeline (Ivančić et al., 2022) to adapt it to our needs (```insert_seq_extract_umi_modif.py```). 
Briefly this script is based on the specific structure of the UMIs integrated in fixed linker sequences. 

The UMI extraction can be performed using the ```extract_umi.sh``` script.

**Usage**
```sh
./extract_umi.sh <DATA_SAM> <DIR_UMI> <nb_barcodes> <adapter_length> <max_error>

Options:
  <DATA_SAM>         | Path to directory of the SAM files obtained after the LTR mapping (barcode${i}_mapping_${a}_SUP.sam)
  <DIR_UMI>          | Path to UMI output directory
  <nb_barcodes>      | Number of barcodes to analyze
  <adapter_length>   | Linker size 
  <max_error>        | Max pattern distance for UMI
```

Outputs in DIR_UMI:
- ```barcode${i}_${LTR}_SUP_fwd.fasta``` : contains the reads in fwd orientation compared to ref LTR sequences
- ```barcode${i}_${LTR}_SUP_rev.fasta``` : contains reads in rev orientation compared to ref LTR sequences
- ```barcode${i}_LTR${a}_UMI.fasta``` : read sequences with identified UMi sequences in the read names
  
## 5- Mapping

## 6- Integration sites extraction

___
## Correspondance

___
## Citation

___
