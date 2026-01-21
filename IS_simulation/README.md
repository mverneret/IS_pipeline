To test the pipeline (basecalling and demultiplexing steps are not included), test files  are available in the ```test_file/``` folder. 
These test files were generated from simulated reads obtained using the scripts available in the ```IS_simulation/scripts/``` folder using domestic goat genome (GCF_001704415.2) and ENTV-2 (3824) as reference sequences. 

**Add read simulation principle**

- ```seq_LTR_linker_CH.fasta``` including the sequences of the LTR5 and 3 (including primers) and the linker sequences (only if you want to simulate reads)

##### Simulate random IS and reads 
```sh
# Simulate 1000 IS and reads with 260bp of flanking host genome taking LTRs from 3824-ENTV2
bash ${WORKDIR}/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r ${WORKDIR}/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l ${WORKDIR}/ref/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o ${WORKDIR}/test_IS/sim \
  > ${WORKDIR}/test_IS/random_IS_CH_1k.log 2>&1
```

