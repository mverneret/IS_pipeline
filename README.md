# Integration Sites pipeline
Pipeline inspired from PCIP-seq ([Artesi et al., 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02307-0)) and INSERT-seq ([Ivančić et al., 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02778-9)) pipelines to identify integration sites (IS) by long read sequencing using Nanopore technology. To do this, the extracted DNA is first fragmented by sonication and the junction between HOST-TARGET sequences (here virus sequences) are amplified by two successive PCRs using specific primers. The resulting reads are filtered to keep only the one including viral and host sequences. The main steps and expected structure of the reads are shown bellow for LTR5. Same principle is applied for LTR3.

<img src="image/Image2.png" width="50%">

___
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
- R (> 4.3.1) **Add R packages**

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
dorado demux --kit-name SQK-NBD114-96 $OUT_BAM_DIR/RUN_${i}_SUP.bam --output-dir $OUTPUT_DIR/demux_both_end --barcode-both-ends --emit-fastq
```
**Add info on the structure/names of the folder and files**

## 3- Filtering steps
Reads were then filtered according to different criteria including :
- Quality > Q20
- Mapping on the LTR5 end or LTR3 start
- Not mapping on provirus without LTR sequences
- Host genome size > 50bp

All these steps can be performed using the ```filter_reads.sh``` script that will deal with LTR5 and LTR3 reads.

**Usage**
```sh
./filter_reads.sh -r <string> -o <string> -f <string> -i <string> -n <string/int> -q <int> -g <int> -a <int> -l <int> -m <int>

Options:
  -r <REF_DIR>         | Path to LTR reference sequences directory
  -o <OUT_DIR>         | Path to output directory
  -f <FASTQ_FILE>      | Path to demultiplexed FASTQ files directory 
  -n <REF_NAME>        | Virus reference name (used for fasta ref files and index names)
  -i <SAMPLE_NAME>     | Sample name (used in prefix for output files)
  -q <MIN_QUALITY>     | Minimum quality for Nanofilt
  -g <LENGTH_GENOME>   | Minimum genome length in reads
  -a <LENGTH_ADAPTOR>  | Adaptor length
  -l <LENGTH_LTR5>     | LTR5 length in reads (with primer)
  -m <LENGHT_LTR3>     | LTR3 length in reads (with primer)
```
Input files in **X folder**:
- Fastq files from demultiplexing
- Fasta sequences of the reference virus without primers (LTR5:"endU3RU5"; LTR3:"stratU3"; INT:"provirus_wo_LTR")
  
Outputs in ```OUT_DIR``` (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- **```barcode${i}_mapping_${a}_SUP.sam``` : Mapping results on start LTR5 or end LTR3**
- ```barcode${i}_mapped_${a}_SUP.sam|fastq``` : Only mapped reads on start LTR5 or end LTR3
- ```barcode${i}_mapped_${a}_mapping_${ref_name}_provirus_wo_LTR_SUP.sam``` : Mapping results on provirus w/o LTR sequences
- ```barcode${i}_${a}_SUP.sam``` : Only non-mapped reads on provirus w/o LTR sequences
- **```barcode${i}_${a}_filtered_size_SUP.sam|fastq``` : Reads after all the steps + size filtering**

## 4- Extract UMI 
In order to remove PCR duplicates for clonality quantification in the Step 6- it was necessary to extract UMI sequences from all the reads. 
We thus modified a python script from the INSERT-seq pipeline (Ivančić et al., 2022) to adapt it to our needs (```insert_seq_extract_umi_modif.py```). 
Briefly this script is based on the specific structure of the UMIs integrated in fixed linker sequences. **Add more details on the extraction of the UMI and on the UMI structure**

The UMI extraction can be performed using the ```extract_umi.sh``` script.

**Usage**
```sh
./extract_umi.sh -s <string> -o <string> -a <int> -e <int> -i <string>

Options:
  -s <DATA_SAM>         | Path to directory of the SAM files obtained after the LTR mapping (barcode${i}_mapping_${a}_SUP.sam)
  -o <DIR_UMI>          | Path to UMI output directory
  -a <ADAPTER_LENGTH>   | Linker size 
  -e <MAX_ERROR>        | Max pattern distance for UMI (allowed error)
  -i <SAMPLE_PREFIX>    | Sample name used as prefix for output files
```
Input files in **X folder**:
- sam files after bowtie2 mapping for LTR5 and 3
  
Outputs in ```DIR_UMI``` (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```barcode${i}_${a}_SUP_fwd.fasta``` : contains the reads in fwd orientation compared to ref LTR sequences
- ```barcode${i}_${a}_SUP_rev.fasta``` : contains reads in rev orientation compared to ref LTR sequences
- ```barcode${i}_LTR${a}_UMI.fasta``` : read sequences with identified UMi sequences in the read names
  
## 5- Mapping
The filtered reads were then mapped on the reference genome concatenated with start LTR5 and end LTR3 sequences in order to detect the HOST-TARGET junctions.
The steps includes :
- Mask the reference genome with the virus and the apparented ERV sequences
- Create the hybrid reference concatenating host genome and LTR sequences
- Map the filtered reads on the hybrid genome using minimap2

The mapping can be performed using the ```mapping.sh``` script.

**Usage**
```sh ./mapping.sh -r REF_DIR -f REF_NAME -t TE_ANNOT -o OUT_DIR -q FASTQ_DIR -c PREFIX_CHR -n VIRUS_NAME -i SAMPLE_PREFIX"

**Add an option to keep/or not the scaffold**
Options:
  -r <REF_DIR>         | Path to directory of the reference files
  -f <REF_NAME>        | Name of the reference genome file in .fasta
  -t <TE_ANNOT>        | Name of the ERV annotation file in .bed
  -o <OUT_DIR>         | Path to directory for output files
  -q <FASTQ_DIR>       | Path to directory containining the input .fastq files
  -c <PREFIX_CHR>      | Prefix of the chromosome names (to keep only chr in the reference and delete scaffolds)
  -n <VIRUS_NAME>      | Name of the virus sequence used in reference for the target
  -i <SAMPLE_PREFIX>   | Sample prefix for output files
```
Input files in **X folder**:
- reference genome in fasta
- ERV annotation in bed
- Virus reference sequences in fasta (with primers)
- fastq files after bowtie2 mapping + filtering

Outputs in ```REF_DIR```  (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```${REF%.fa}_noscaffold_masked.fa``` : Masked reference genome without scaffolds
- ```${REF%.fa}_masked_${a}_withprimer.fa``` : Masked reference genome without scaffolds + LTR sequence 
- ```${REF%.fa}_masked_${a}_withprimer.mmi``` : Indexed hybrid reference

Outputs in ```OUT_DIR```:
- ```barcode${i}_LTR${a}_mapped_${REF}_SUP.sam|bam|paf``` : minimap2 output file
- ```barcode${i}_LTR${a}_mapped_${REF}_sorted_SUP.bam```

## 6- Integration sites extraction
After the mapping, the goal is identify the different integration sites using reads at the junction between LTR sequences and host genome. The main steps are:
- Keep reads that mapped on the host genome and on the LTR sequences
- Get the Integration Sites (IS) corresponding to the HOST-LTR junction and ShearSites (ShS) corresponding to HOST-LINKER junction
- Create ShS groups clustering reads with ShS < maxgapShS + UMI group using ```UMI_clustering_hamming_ref.py``` if reads have same UMI (+/- x mismatches)
- Remove PCR duplicates according to the IS, ShS and UMI groups
- Merge LTR5 and LTR3 information 
- Quantify clonal abundancy (for each IS nb of sister cells without PCR duplicates)

All these steps can be performed using ```run_IS.sh``` script that will call different R functions.

**Usage**
```sh
./run_IS.sh <sample_name> <R_package_path> <out_path> <input_paf_path> <input_UMI_path> <assembly> <targetName_LTR5> \
    <lengthTarget_LTR5> <targetName_LTR3> <lengthTarget_LTR3> <maxgapIS> <mapq_val> <nb_read_per_umi> <maxgapShS> <mms> <threshold_raw> <win_merge>

Options:
  <sample_name>         | Name of the sample (ex: barcode01)
  <R_package_path>      | Path to directory of the different R functions (ex: ~/script/Rpackages/)
  <out_path>            | Path to directory for output files (ex: ~/results/R_clonality/)
  <input_paf_path>      | Path to directory of the paf input files obtained after the step 5 (ex: ~/results/mapping/paf/)
  <input_UMI_path>      | Path to directory of the UMI input files obtained after the step 4 (ex: ~/results/extract_UMI/)
  <assembly>            | Name of the reference assembly
  <targetName_LTR5>     | Name of the LTR5 virus sequence (used for the mapping)
  <lengthTarget_LTR5>   | Total length of the LTR5 sequence
  <targetName_LTR3>     | Name of the LTR3 virus sequence (used for the mapping)
  <lengthTarget_LTR3>   | Total length of the LTR3 sequence
  <maxgapIS>            | Maximal distance between IS to merge
  <mapq_val>            | Minimum quality of mapping
  <nb_read_per_umi>     | Minimum number of reads to conserve a UMI group
  <maxgapShS>           | Maximal distance between ShS to merge
  <mms>                 | Maximum number of mismatchs to regroup UMI in a group
  <threshold_raw>       | Minimum number of raw reads by IS to keep it
  <win_merge>           | Maximal distance between LTR5 and LTR3 IS to merge
```

Outputs in ```out_path``` (i = barcode number; a= LTR5 (startU3) or LTR3 (endU3RU5)):
- ```*merged*``` : Each reads with all their info + attributed UMI group
- ```*positionreads*``` : Each reads with all their info + UMI, ShS and IS group
- ```*countedreadsLTR3|LTR5*``` : Grouped IS with coordinates and nb of reads
- ```*countedreadsLTR5LTR3*``` : Grouped IS with coordinates and nb of reads after LTR5 and LTR5 merge
- ```*clonalityResults*.txt``` : Final IS results with clonality %

**Add info about test dataset and scripts + parameters file**
___
## Correspondance
Benjamin Riocreux-Verney,
Jocelyn Turpin (jocelyn.turpin@univ-lyon1.fr)

## Citation
___
