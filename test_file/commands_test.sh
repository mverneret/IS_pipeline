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

####### Mask and filter out scaffolds from the host reference genome #######
bedtools maskfasta -fi GCF_001704415.2.fa -bed CH_annotation_II-5_ARS1.2.bed -fo GCF_001704415.2_masked.fa
prefix_chr="NC"
awk -v prefix_chr="^>${prefix_chr}" ' BEGIN { keep=0 } /^>/ { keep = ($0 ~ prefix_chr) } keep { print } ' GCF_001704415.2_masked.fa > ARS12_noscaffold_masked.fa

####### Simulate IS and reads #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/data/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/sim \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k.log 2>&1

####### Align the random reads on LTR sequences and separate LTR5 and LTR3 #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/filter_reads.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/bowtie2 \
  -f /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/sim/random_IS_CH_1k_reads.fq \
  -i random_IS_CH_1k \
  -n 3824 \
  -q 20 \
  -g 50 \
  -a 57 \
  -l 63 \
  -m 188 \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_bowtie2.log 2>&1

####### Extract UMIs from the mapped reads #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/extract_umi.sh \
  -s /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/bowtie2 \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/extract_UMI \
  -a 57 \
  -e 0 \
  -i random_IS_CH_1k \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_extract_umi.log 2>&1

####### Mapping of the reads on the virus + host ref genomes #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/mapping.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref \
  -f ARS12_noscaffold_masked.fa \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/mapping \
  -q /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/bowtie2 \
  -c NC \
  -n 3824 \
  -i random_IS_CH_1k \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_mapping.log 2>&1

####### Identify integration sites #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/R/run_IS.sh \
  -i random_IS_CH_1k \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/R/Rpackage/ \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/R_clonality/ \
  -p /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/mapping/paf/ \
  -u /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/extract_UMI/ \
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
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_IS.log 2>&1
