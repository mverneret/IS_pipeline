####### VM and tools install #######
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


git clone https://github.com/mverneret/IS_pipeline.git
WORKDIR="$(pwd)" 
##Ex on VM : "/home/ubuntu/data/mydatalocal"

mkdir ${WORKDIR}/ref

#### What files in ref/ and how download the reference genome
cd ref
####### Mask and filter out scaffolds from the host reference genome #######
bedtools maskfasta -fi GCF_001704415.2.fa -bed CH_annotation_II-5_ARS1.2.bed -fo GCF_001704415.2_masked.fa
prefix_chr="NC"
awk -v prefix_chr="^>${prefix_chr}" ' BEGIN { keep=0 } /^>/ { keep = ($0 ~ prefix_chr) } keep { print } ' GCF_001704415.2_masked.fa > ARS12_noscaffold_masked.fa

mkdir ${WORKDIR}/test_IS
cd test_IS

mkdir R_clonality bowtie2 extract_UMI mapping sim
