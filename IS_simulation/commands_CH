bash /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/scripts/simulate_reads.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref/ARS12_noscaffold_masked.fa \
  -n 1000 \
  -f 260 \
  -l /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/data/seq_LTR_linker_CH.fasta \
  -d random_IS_CH_1k \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/sim \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k.log 2>&1


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


bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/extract_umi.sh \
  -s /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/bowtie2 \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/extract_UMI \
  -a 57 \
  -e 0 \
  -i random_IS_CH_1k \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_extract_umi.log 2>&1

bash /home/ubuntu/data/mydatalocal/IS_pipeline/scripts/mapping.sh \
  -r /home/ubuntu/data/mydatalocal/IS_pipeline/test/ref \
  -f ARS12_noscaffold_masked.fa \
  -o /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/mapping \
  -q /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/bowtie2 \
  -c NC \
  -n 3824 \
  -i random_IS_CH_1k \
  > /home/ubuntu/data/mydatalocal/IS_pipeline/IS_simulation/random_IS_CH_1k_mapping.log 2>&1


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
