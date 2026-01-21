# IS and reads simulation
To test the pipeline (basecalling and demultiplexing steps are not included), test files  are available in the ```test_file/``` folder. 
These test files were generated from simulated reads obtained using the scripts available in the ```IS_simulation/scripts/``` folder using domestic goat genome (GCF_001704415.2) and ENTV-2 (3824) as reference sequences. 

Reads are simulated using several steps:
- Generate random IS sites taking the chromosome size into account using ```random_sites.py```
- Extract flanking region of the random IS : either in 5' or in 3' randomly using ```generate_flank_bed.py```
- Build simulated reads simulating UMIs and using linker, LTR and host flanking genome sequences with ```final_random_reads.py```
- Convert fasta to fastq file

The structure of the reads is build as bellow:


All the simulation can be excecuted using the ```simulate_reads.sh``` script.

### Usage
```sh
./simulate_reads.sh -r <string> -n <int> -f <int> -l <string> -d <string> -o <string>

Options:
  -r <REF_FASTA>         | Path to reference genome after masking and removing scaffolds
  -n <NUM_SITES>         | Number of integration sites to generate
  -f <FLANK_SIZE>        | Fixed flanking genome size (bp) 
  -l <LTR_LINKER_FASTA>  | Path to fasta file containing the linker + LTR sequences
  -d <OUT_PREFIX>        | Output prefix
  -o <OUT_DIR>           | Path to outputs
```

Inputs files:
- ```ref_genome_noscaffold_masked.fa``` : prepared host reference genome 
- ```seq_LTR_linker.fasta``` : the linker sequences and the virus LTR5 and LTR3 sequences (including primers)

Output files:
- ```${PREFIX}.bed``` : position of the simulated IS in bed format
- ```${PREFIX}_reads.fq``` : fastq files of the simulated reads

### Example of run to produce the test dataset
```sh
WORKDIR="$(pwd)" 
##Ex on VM : "/home/ubuntu/data/mydatalocal"

# Domestic goat is used as reference genome (GCF_001704415.2). Check the 3- step of the test_file/README.md for genome download and preparation
# LTRs sequences from 3824-ENTV2 are used and present in the seq_LTR_linker.fasta file for simulation
# Simulate 1000 IS and reads with 260bp of flanking goat genome
bash ${WORKDIR}/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r ${WORKDIR}/IS_pipeline/test_file/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l ${WORKDIR}/IS_pipeline/test_file/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o ${WORKDIR}/IS_pipeline/test_file \
  > ${WORKDIR}/IS_pipeline/test_file/random_IS_CH_1k.log 2>&1
```
