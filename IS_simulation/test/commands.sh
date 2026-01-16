####### Simulate IS and reads #######
bash /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/data/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/sim \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k.log 2>&1
