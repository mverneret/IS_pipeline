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
  -n <ref_name>        | Virus reference name (used for fasta ref files and index names) 
  -q <min_quality>     | Minimum quality for Nanofilt
  -g <length_genome>   | Minimum genome length in reads
  -a <length_adaptor>  | Adaptor length
  -l <length_LTR5>     | LTR5 length in reads (with primer)
  -m <length LTR3>     | LTR3 length in reads (with primer)
  -b <nb_barcode>      | Number of barcodes to analyze
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

Outputs in DIR_UMI (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```barcode${i}_${a}_SUP_fwd.fasta``` : contains the reads in fwd orientation compared to ref LTR sequences
- ```barcode${i}_${a}_SUP_rev.fasta``` : contains reads in rev orientation compared to ref LTR sequences
- ```barcode${i}_LTR${a}_UMI.fasta``` : read sequences with identified UMi sequences in the read names
  
## 5- Mapping
The filtered reads were then mapped on the reference genome concatenated with start LTR5 and end LTR3 sequences in order to detect the HOST-TARGET junctions.
The first step is then to mask the reference genome with the virus and the apparented ERV sequences, then create the hybrid reference and map the filtered reads on it using minimap2.

The mapping can be performed using the ```mapping.sh``` script.

**Usage**
```sh
./mapping.sh <REF_DIR> <REF_name> <TE_annot> <OUT_DIR> <FASTQ_DIR> <prefix_chr> <virus_name> <nb_barcodes>

Options:
  <REF_DIR>         | Path to directory of the reference genome
  <REF_name>          | Name of the reference genome file in .fasta
  <TE_annot>      | Name of the ERV annotation file in .bed
  <OUT_DIR>   | Path to directory for output file
  <FASTQ_DIR>        | Path to directory contaiining the input .fastq files
  <prefix_chr>   | Préfix of the chromosome names (to keep only chr in the refernece and delete scaffolds)
  <virus_name>  | Name of the virus/sequence used in reference for the target
  <nb_barcodes>  | Number of barcodes to analyze
```

Outputs in REF_DIR (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```${REF%.fa}_noscaffold_masked.fa``` : Masked reference genome without scaffolds
- ```${REF%.fa}_masked_${a}_withprimer.fa``` : Masked reference genome without scaffolds + LTR sequence 
- ```${REF%.fa}_masked_${a}_withprimer.mmi``` : Indexed hybrid reference

Outputs in OUT_DIR:
- ```barcode${i}_LTR${a}_mapped_${REF}_SUP.sam|bam``` : minimap2 output file
- ```barcode${i}_LTR${a}_mapped_${REF}_sorted_SUP.bam```

## 6- Integration sites extraction

___
## Correspondance

___
## Citation

___
