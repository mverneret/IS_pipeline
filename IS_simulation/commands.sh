############################ Prepare the reference genome #####################################
#mask genome with position of ERV II.5 on ARS1.2
bedtools maskfasta -fi GCF_001704415.2.fa -bed CH_annotation_II-5_ARS1.2.bed -fo GCF_001704415.2_masked.fa
#remove scafold
awk '
BEGIN {header=""; sequence=""}
/^>/ {
#rajouter || header ~ /^>ENTV/
    if (header != "" && (header ~ /^>NC/) {
        print header
        print sequence
    }
    header=$0
    sequence=""
    next
}
{
    sequence = sequence $0 "\n"
}
END {
#rajouter || header ~ /^>ENTV/
    if (header ~ /^>NC/) {
        print header
        print sequence
    }
}
' GCF_001704415.2_masked.fa > ${DATA_DIR}/ARS12_noscaffold_masked.fa

############################## Produce the simulated reads #####################################

#produce a bed file with random integration sites taking the chromosome size into consideration (wieghted chr depending on the size)
python3 random_sites.py --fasta ARS12_noscaffold_masked.fa --output random_int_1M.bed --num-sites 1000000

sort -k1,1 -2,2n random_int_1M.bed > random_int_1M_sorted.bed

#number of integration sites by chr
echo "Number of integration sites per chromosome:"
cut -f1 random_int_1M_sorted.bed | sort | uniq -c | sort -k2

#extract randomly the 5prime of 3prime flanking sequences of the integration sites (260bp)
python3 generate_flank_bed.py random_int_1M_sorted.bed random_int_1M_sorted_flank.bed 260

bedtools getfasta -fi ARS12_noscaffold_masked.fa -bed random_int_1M_sorted_flank.bed -fo random_int_1M_sorted_flank.fasta -nameOnly

#concatenate linker, LTR and flanking sequences to produce the simulated reads
python3 build_final_random_reads.py random_int_1M_sorted_flank.fasta seq_LTR_linker.fasta random_int_1M_sorted_flank_reads.fasta

#count number of reads for each combination
echo "Number reads for 5prime_sens_+:"
grep "5prime_sens_+" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 5prime_antisens_+:"
grep "5prime_antisens_+" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 5prime_sens_-:"
grep "5prime_sens_-" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 5prime_antisens_-:"
grep "5prime_antisens_-" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 3prime_sens_+:"
grep "3prime_sens_+" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 3prime_antisens_+:"
grep "3prime_antisens_+" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 3prime_sens_-:"
grep "3prime_sens_-" random_int_1M_sorted_flank_reads.fasta  | wc -l

echo "Number reads for 3prime_antisens_-:"
grep "3prime_antisens_-" random_int_1M_sorted_flank_reads.fasta  | wc -l


################# Align the random reads on LTR sequences and separate LTR5 and LTR3 #####################################"
cd botwie2/

#separate LTR5 and LTR3
awk '/^>/{keep=0} /^>.*(5prime_sens|3prime_antisens)/{keep=1} keep' random_int_1M_sorted_flank_reads.fasta > random_LTR5_SUP.fasta
awk '/^>/{keep=0} /^>.*(3prime_sens|5prime_antisens)/{keep=1} keep' random_int_1M_sorted_flank_reads.fasta > random_LTR3_SUP.fasta

#convert fasta in fastq
seqtk seq -F 'I' random_LTR3_SUP.fasta > random_LTR3_SUP.fq
seqtk seq -F 'I' random_LTR5_SUP.fasta > random_LTR5_SUP.fq

#index the LTR references
/beegfs/data/soft_legacy/bowtie2-2.1.0/bowtie2-build ENTV2_startU3_withprimer.fa ENTV2_startU3_index
/beegfs/data/soft_legacy/bowtie2-2.1.0/bowtie2-build ENTV2_endU3RU5_withprimer.fa ENTV2_endU3RU5_index

#map reads on LTR references
/beegfs/data/soft_legacy/bowtie2-2.1.0/bowtie2 --sensitive-local -x ENTV2_endU3RU5_index -U random_LTR3_SUP.fq -S random_LTR3_SUP.sam -N 1
/beegfs/data/soft_legacy/bowtie2-2.1.0/bowtie2 --sensitive-local -x ENTV2_startU3_index -U random_LTR5_SUP.fq -S random_LTR5_SUP.sam -N 1

#convert sam in fastq
samtools fastq random_LTR3_SUP.sam > random_LTR3_SUP.fastq
samtools fastq random_LTR5_SUP.sam > random_LTR5_SUP.fastq

#extract fwd and rev mapped reads
samtools view -F 20 random_LTR3_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > random_LTR3_SUP_fwd.fasta
samtools view -f 16 -F 4 random_LTR3_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > random_LTR3_SUP_rev.fasta

samtools view -F 20 random_LTR5_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > random_LTR5_SUP_fwd.fasta
samtools view -f 16 -F 4 random_LTR5_SUP.sam | grep -v '^@' | awk '{print ">" $1 "\n" $10}' > random_LTR5_SUP_rev.fasta

########################## Extract UMIs from the mapped reads ##########################################
cd extract_UMI/

python3 ../../results/script/insert_seq_extract_umi_modif3.py --max-error 0 --adapter-length 57 -o random_LTR5_UMI.fasta -ifwd random_LTR5_SUP_fwd.fasta -irev random_LTR5_SUP_rev.fasta
python3 ../../results/script/insert_seq_extract_umi_modif3.py --max-error 0 --adapter-length 57 -o random_LTR3_UMI.fasta -ifwd random_LTR3_SUP_fwd.fasta -irev random_LTR3_SUP_rev.fasta


######################### Mapping of the random reads on the goat ref genome ##########################
cd mapping/

#concatenate goat ref genome with LTR sequences
cat ../ARS12_noscaffold_masked.fa ENTV2_endU3RU5_withprimer.fa > ARS12_noscaffold_masked_U5_withprimers.fa
cat ../ARS12_noscaffold_masked.fa ENTV2_startU3_withprimer.fa > ARS12_noscaffold_masked_U3_withprimers.fa

#index the new hybrid ref
minimap2 -d ARS12_noscaffold_masked_U5_withprimer.mmi ARS12_noscaffold_masked_U5_withprimer.fa
minimap2 -d ARS12_noscaffold_masked_U3_withprimer.mmi ARS12_noscaffold_masked_U3_withprimer.fa

#map the random reads previously mapped on LTR on the goat and LTR reference
minimap2 -ax map-ont -t 8 ARS12_noscaffold_masked_U5_withprimer.mmi random_LTR3_SUP.fastq > random_LTR3_mapped_ARS12_SUP.sam
minimap2 -ax map-ont -t 8 ARS12_noscaffold_masked_U3_withprimer.mmi random_LTR5_SUP.fastq > random_LTR5_mapped_ARS12_SUP.sam

#convert sam to paf
paftools.js sam2paf random_LTR3_mapped_ARS12_SUP.sam > random_LTR3_mapped_ARS12_SUP.paf
paftools.js sam2paf random_LTR5_mapped_ARS12_SUP.sam > random_LTR5_mapped_ARS12_SUP.paf

#convert sam to bam
samtools view -@ 4 -Sb random_LTR5_mapped_ARS12_SUP.sam > random_LTR5_mapped_ARS12_SUP.bam
samtools view -@ 4 -Sb random_LTR3_mapped_ARS12_SUP.sam > random_LTR3_mapped_ARS12_SUP.bam
#index the bam
samtools sort -@ 4 random_LTR5_mapped_ARS12_SUP.bam > random_LTR5_mapped_ARS12_sorted_SUP.bam
samtools sort -@ 4 random_LTR3_mapped_ARS12_SUP.bam > random_LTR3_mapped_ARS12_sorted_SUP.bam
samtools index -@ 4 random_LTR3_mapped_ARS12_sorted_SUP.bam
samtools index -@ 4 random_LTR5_mapped_ARS12_sorted_SUP.bam
