# IS and reads simulation
To test the pipeline (basecalling and demultiplexing steps are not included), test files  are available in the ```test_file/``` folder. 
These test files were generated from simulated reads obtained using the scripts available in the ```IS_simulation/scripts/``` folder using domestic goat genome (GCF_001704415.2) and ENTV-2 (3824) as reference sequences. 

**Add read simulation principle**


### Run the script
```sh
WORKDIR="$(pwd)" 
##Ex on VM : "/home/ubuntu/data/mydatalocal"

# Simulate 1000 IS and reads with 260bp of flanking goat genome and taking LTRs from 3824-ENTV2
bash ${WORKDIR}/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r ${WORKDIR}/IS_pipeline/test_file/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l ${WORKDIR}/IS_pipeline/test_file/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o ${WORKDIR}/IS_pipeline/test_file \
  > ${WORKDIR}/IS_pipeline/test_file/random_IS_CH_1k.log 2>&1
```

Inputs files:
- ```ARS12_noscaffold_masked.fa``` : prepared domestic goat reference genome (GCF_001704415.2). Check the 3- step of the ```test_file/README.md``` for genome download and preparation
- ```seq_LTR_linker_CH.fasta``` : the linker sequences and the virus (here 3824-ENTV2) LTR5 and LTR3 sequences (including primers)

Output files:
- ```random_IS_CH_1k.bed``` : position of the simulated IS in bed format
- ```random_IS_CH_1k_reads.fq``` : fastq files of the simulated reads
