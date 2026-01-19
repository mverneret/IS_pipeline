#!/bin/bash

####### Simulate IS and reads #######
bash ${WORKDIR}/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r ${WORKDIR}/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l ${WORKDIR}/ref/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o ${WORKDIR}/test_IS/sim \
  > ${WORKDIR}/test_IS/random_IS_CH_1k.log 2>&1
  
####### Align the random reads on LTR sequences and separate LTR5 and LTR3 #######
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

####### Extract UMIs from the mapped reads #######
bash ${WORKDIR}/IS_pipeline/scripts/extract_umi.sh \
  -s ${WORKDIR}/test_IS/bowtie2 \
  -o ${WORKDIR}/test_IS/extract_UMI \
  -a 57 \
  -e 0 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_extract_umi.log 2>&1

####### Mapping of the reads on the virus + host ref genomes #######
bash ${WORKDIR}/IS_pipeline/scripts/mapping.sh \
  -r ${WORKDIR}/ref \
  -f ARS12_noscaffold_masked.fa \
  -o ${WORKDIR}/test_IS/mapping \
  -q ${WORKDIR}/test_IS/bowtie2 \
  -n 3824 \
  -i random_IS_CH_1k \
  > ${WORKDIR}/test_IS/random_IS_CH_1k_mapping.log 2>&1

####### Identify integration sites #######
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
