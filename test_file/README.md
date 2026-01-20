# Test dataset and reads simulation
To test the pipeline (basecalling and demultiplexing steps are not included), test files  are available in the ```test_file/``` folder. 
These test files were generated from simulated reads obtained using the scripts available in the ```IS_simulation/scripts/``` folder using domestic goat genome (GCF_001704415.2) and ENTV-2 (3824) as reference sequences. 

**Add read simulation principle**

## Running tutorial
#### 1- VM and tools install
```sh
sudo apt update
sudo apt install r-base r-base-dev
sudo apt install python3-edlib python3-pysam python3-tqdm python3-pandas python3-biopython
sudo apt install samtools bedtools seqtk bowtie2 nanofilt
#install minimap2 and paftools.js
git clone https://github.com/lh3/minimap2
cd minimap2 && make
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8  
#install R packages
sudo apt install libbz2-dev liblzma-dev zlib1g-dev
#in R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
install.packages("changepoint")
install.packages("tidyverse")
BiocManager::install("GenomicAlignments")
BiocManager::install("stringdist")
```

#### 2- Clone the IS_pipeline repository
```sh
git clone https://github.com/mverneret/IS_pipeline.git
```

#### 3- Create the ref/ folder and reference files
```sh
WORKDIR="$(pwd)" 
##Ex on VM : "/home/ubuntu/data/mydatalocal"

mkdir ${WORKDIR}/ref
cd ref
```

This ```ref/``` folder must contain the host reference genome and the virus reference sequences (available in the ```test_file/ref``` folder of this repo) :
- ```ARS12_reference_genome.fa```
- ```3824_startU3.fa``` + ```3824_startU3_withprimers.fa``` : LTR5 virus reference sequences
- ```3824_endU3RU5.fa``` + ```3824_endU3RU5_withprimers.fa``` : LTR3 virus reference sequences
- ```3824_provirus_wo_LTR.fa``` : INT part without LTR virus reference sequences
- ```seq_LTR_linker_CH.fasta``` including the sequences of the LTR5 and 3 (including primers) and the linker sequences (only if you want to simulate reads)

The reference host genome can be downloaded as bellow :
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.2_ARS1.2/GCF_001704415.2_ARS1.2_genomic.fna.gz
gzip -d GCF_001704415.2_ARS1.2_genomic.fna.gz
####### Mask and filter out scaffolds from the host reference genome #######
bedtools maskfasta -fi GCF_001704415.2_ARS1.2_genomic.fna -bed Annotation_ERV_ARS1.2_GCF_001704415.2_family_II-5.bed -fo GCF_001704415.2_masked.fa
prefix_chr="NC"
awk -v prefix_chr="^>${prefix_chr}" ' BEGIN { keep=0 } /^>/ { keep = ($0 ~ prefix_chr) } keep { print } ' GCF_001704415.2_masked.fa > ARS12_noscaffold_masked.fa
```

#### 4- Create de results folder
Different folders must be created for the outputs of the pipeline
- ```sim/``` : results of the IS and reads simulation
- ```bowtie2/``` : results of the 3-Filtering steps
- ```extract_umi/```: results of the 4-Extract UMI
- ```mapping/``` : results of the 5-Mapping
- ```Rclonality/``` : results of the 6-Integration sites extraction
    
```sh
mkdir ${WORKDIR}/test_IS
cd test_IS
mkdir R_clonality bowtie2 extract_UMI mapping sim
```

#### 5- Test if the pipeline is working

####### Simulate IS and reads #######
```sh
bash ${WORKDIR}/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r ${WORKDIR}/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l ${WORKDIR}/ref/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o ${WORKDIR}/test_IS/sim \
  > ${WORKDIR}/test_IS/random_IS_CH_1k.log 2>&1
```

####### Align the random reads on LTR sequences and separate LTR5 and LTR3 #######
```sh
bash ${WORKDIR}/IS_pipeline/scripts/filter_reads.sh \
  -r ${WORKDIR}/ref \
  -o ${WORKDIR}/test_IS/bowtie2 \
  -f ${WORKDIR}/test_IS/sim/random_IS_CH_1k_reads.fq \
  -i random_IS_CH_1k \
  -n 3824 \
  -q 20 \
  -g 50 \
  -a 57 \
  -l 63 \
  -m 188 \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_bowtie2.log 2>&1
```

####### Extract UMIs from the mapped reads #######
```sh
bash ${WORKDIR}/IS_pipeline/scripts/extract_umi.sh \
  -s ${WORKDIR}/test_IS/bowtie2 \
  -o ${WORKDIR}/test_IS/extract_UMI \
  -a 57 \
  -e 0 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_extract_umi.log 2>&1
```

####### Mapping of the reads on the virus + host ref genomes #######
```sh
bash ${WORKDIR}/IS_pipeline/scripts/mapping.sh \
  -r ${WORKDIR}/ref \
  -f ARS12_noscaffold_masked.fa \
  -o ${WORKDIR}/test_IS/mapping \
  -q ${WORKDIR}/test_IS/bowtie2 \
  -n 3824 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_mapping.log 2>&1
```

####### Identify integration sites #######
```sh
bash ${WORKDIR}/IS_pipeline/scripts/R/run_IS.sh \
  -i random_IS_CH_1k \
  -r ${WORKDIR}/IS_pipeline/scripts/R/Rpackage/ \
  -o ${WORKDIR}/test_IS/R_clonality/ \
  -p ${WORKDIR}/test_IS/mapping/paf/ \
  -u ${WORKDIR}/test_IS/extract_UMI/ \
  -a "ARS12_noscaffold_masked.fa" \
  -5 "3824_startU3" \
  -L 63 \
  -3 "3824_endU3RU5" \
  -l 188 \
  -g 15 \
  -q 20 \
  -n 1 \
  -s 10 \
  -m 2 \
  -t 1 \
  -w 25 \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_IS.log 2>&1
```
