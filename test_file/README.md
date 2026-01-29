# Install and test dataset
The pipeline can be tested using the test dataset composed of simulated reads obtained using ARS1.2 domestic goat genome (GCF_001704415.2) and ENTV-2 (strain FR3824, PP669281.1) as reference sequences.\
Check on ```IS_simulation/README.rmd``` for more information on how the simulated reads were produced.

All the steps to install and set up the pipeline are detailed bellow. 

### 1- Clone the IS_pipeline repository
```sh
git clone https://github.com/mverneret/IS_pipeline.git
```

### 2- Install dependancies

**Using mamba**

- Re-create the mamba environment
```sh
mamba env create -f IS_pipeline.yml
```

- **OR** create the mamba environment directly
```sh
mamba create -n IS_pipeline \
  python=3.10.14 \
  r-base=4.3.3 \
  bowtie2=2.5.2 \
  nanofilt=2.8.0 \
  samtools=1.16.1 \
  bedtools=2.30.0 \
  seqtk=1.4 \
  minimap2=2.26 \
  edlib \
  pysam \
  tqdm \
  pandas \
  biopython \
  r-biocmanager \
  r-stringdist \
  r-changepoint \
  r-tidyverse \
  bioconductor-genomicranges \
  bioconductor-rsamtools \
  bioconductor-genomicalignments \
  -c conda-forge -c bioconda
```

- Activate the mamba environment:
```sh
mamba activate IS_pipeline
```

**Using sudo rights**

```sh
sudo apt update
#install R
sudo apt install r-base r-base-dev
#install python3 packages
sudo apt install python3-edlib python3-pysam python3-tqdm python3-pandas python3-biopython
#install other tools needed in the pipeline
sudo apt install samtools bedtools seqtk bowtie2 nanofilt

#install minimap2 and paftools.js
git clone https://github.com/lh3/minimap2
cd minimap2 && make
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8
#replace by the minimap2 path
export PATH="$PATH:/home/ubuntu/minimap2:/home/ubuntu/minimap2/misc"
source ~/.bashrc

#install R packages
# ---------- system dependencies ----------
sudo apt update
sudo apt install -y \
  build-essential \
  pkg-config \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  zlib1g-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libbz2-dev \
  liblzma-dev

# ---------- R packages (CRAN + Bioconductor) ----------
#in R
install.packages('BiocManager', repos='https://cloud.r-project.org')
BiocManager::install(c('GenomicRanges', 'Rsamtools', 'GenomicAlignments'))
install.packages(c('stringdist', 'changepoint'), repos = 'https://cloud.r-project.org')
# ---------- test libraries ----------
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(stringdist)
library(changepoint)

#Install specifically tidyverse
sudo apt update
sudo apt install -y \
  libfontconfig1-dev \
  libfreetype6-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  pkg-config
# In R
install.packages("tidyverse")
library(tidyverse)
```

### 3- The ref/ folder and reference files
```sh
WORKDIR="$(pwd)" 
##Ex on VM : "/home/ubuntu/data/mydatalocal"

cd ${WORKDIR}/IS_pipeline/test_file/ref
```

This ```ref/``` folder must contain the host reference genome and the virus reference sequences :
- ```ARS12_noscaffold_masked.fa``` : Domestic goat ARS1.2 reference genome (GCF_001704415.2)
- ```3824_startU3.fa``` + ```3824_startU3_withprimers.fa``` : ENTV-2 (strain FR3824, PP669281.1) LTR5 reference sequences
- ```3824_endU3RU5.fa``` + ```3824_endU3RU5_withprimers.fa``` : ENTV-2 (strain FR3824, PP669281.1) LTR3 reference sequences
- ```3824_provirus_wo_LTR.fa``` : ENTV-2 (strain FR3824, PP669281.1) provirus without LTR reference sequences

The reference host genome can be downloaded as bellow :
```sh
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.2_ARS1.2/GCF_001704415.2_ARS1.2_genomic.fna.gz
gzip -d GCF_001704415.2_ARS1.2_genomic.fna.gz
####### Mask and filter out scaffolds from the host reference genome #######
bedtools maskfasta -fi GCF_001704415.2_ARS1.2_genomic.fna -bed CH_annotation_II-5_ARS1.2.bed -fo GCF_001704415.2_masked.fa
prefix_chr="NC"
awk -v prefix_chr="^>${prefix_chr}" ' BEGIN { keep=0 } /^>/ { keep = ($0 ~ prefix_chr) } keep { print } ' GCF_001704415.2_masked.fa > ARS12_noscaffold_masked.fa
rm GCF_001704415.2*
```

### 4- Create the results folders
Different folders must be created for the outputs of the pipeline
- ```bowtie2/``` : results of the 3-Filtering steps
- ```extract_umi/```: results of the 4-Extract UMI
- ```mapping/``` : results of the 5-Mapping
- ```Rclonality/``` : results of the 6-Integration sites extraction
    
```sh
mkdir ${WORKDIR}/test_IS
cd ${WORKDIR}/test_IS
mkdir R_clonality bowtie2 extract_UMI mapping
```

### 5- Pipeline test

##### Align the reads on LTR sequences and separate LTR5 and LTR3 
```sh
# Keep reads with quality Q > 20
# Reads aligned on ENTV2-FR3824 LTR5 and LTR3 reference sequences and keep those with
# Read with length > length linker + mapped genome + LTR
# For LTR5 read must have length > 57+50+63
# For LTR3 read must have length > 57+50+188
bash ${WORKDIR}/IS_pipeline/scripts/filter_reads.sh \
  -r ${WORKDIR}/IS_pipeline/test_file/ref \
  -o ${WORKDIR}/test_IS/bowtie2 \
  -f ${WORKDIR}/IS_pipeline/IS_simulation/test_file/random_IS_CH_1k_reads.fq \
  -i random_IS_CH_1k \
  -n 3824 \
  -q 20 \
  -g 50 \
  -a 57 \
  -l 63 \
  -m 188 \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_bowtie2.log 2>&1
```

##### Extract UMIs from the mapped reads
```sh
# UMIs are located in the 57bp linker and 0 mismatchs are allowed between UMI sequences and the expected pattern
bash ${WORKDIR}/IS_pipeline/scripts/extract_umi.sh \
  -s ${WORKDIR}/test_IS/bowtie2 \
  -o ${WORKDIR}/test_IS/extract_UMI \
  -a 57 \
  -e 0 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_extract_umi.log 2>&1
```

##### Mapping of the reads on the virus + host ref genomes
```sh
# Reads are mapped on ARS1.2 (GCF_001704415.2) domestic goat genome concatenated with ENTV2 - FR3824 LTR reference sequences 
bash ${WORKDIR}/IS_pipeline/scripts/mapping.sh \
  -r ${WORKDIR}/IS_pipeline/test_file/ref \
  -f ARS12_noscaffold_masked.fa \
  -o ${WORKDIR}/test_IS/mapping \
  -q ${WORKDIR}/test_IS/bowtie2 \
  -n 3824 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_mapping.log 2>&1
```

##### Identify integration sites
```sh
bash ${WORKDIR}/IS_pipeline/scripts/R/run_IS.sh \
  -i random_IS_CH_1k \
  -r ${WORKDIR}/IS_pipeline/scripts/R/Rpackage/ \
  -o ${WORKDIR}/test_IS/R_clonality/ \
  -p ${WORKDIR}/test_IS/mapping/paf/ \
  -u ${WORKDIR}/test_IS/extract_UMI/ \
  -a "ARS12_noscaffold_masked.fa" \
  -5 "3824_startU3" \
  -3 "3824_endU3RU5" \
  -g 15 \
  -q 20 \
  -n 1 \
  -s 10 \
  -m 2 \
  -t 1 \
  -w 25 \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_IS.log 2>&1
```

All the ouputs of the pipeline on the test dataset are given as an exemple in the ```test_file/results``` folder.
