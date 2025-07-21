#!/bin/bash
#Need to be able to lauch direclty R script via command line

#SBATCH --job-name=R_clonality_barcode01
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH -p workq
#SBATCH --output=run_IS_barcode01.out

module load statistics/R/4.3.1

if [ "$#" -ne 17 ]; then
    echo "Usage: $0 >sample_name> <R_package_path> <out_path> <input_paf_path> <input_UMI_path> <assembly> <targetName_LTR5> n\
    <lengthTarget_LTR5> <targetName_LTR3> <lengthTarget_LTR3> <maxgapIS> <mapq_val> <nb_read_per_umi> <maxgapShS> <mms> <threshold_raw> <win_merge>"
    exit 1
fi

sample_name="$1"
R_package_path="$2"
out_path="$3"
input_paf_path="$4"
input_umi_path="$5"
assembly="$6"

###LTR5###
targetName_LTR5="$7"
lengthTarget_LTR5="$8"

###LTR3###
targetName_LTR3="$9"
lengthTarget_LTR3="$10"

maxgapIS="$11"
mapq_val="$12"
nb_read_per_umi="$13"
maxgapShS="$14"
mms="$15"
threshold_raw="$16"
win_merge="$17"

#Run the first part of the R script : keep chimeric reads, get junctions, add UMI
echo "Starting first R script..."
Rscript /home/mverneret/work/integration/script/run_IS_script_1.R $R_package_path $sample_name $out_path $input_paf_path \
$input_umi_path $targetName_LTR5 $lengthTarget_LTR5 $targetName_LTR3 $lengthTarget_LTR3 $mapq_val $assembly
echo "First R script completed."

#UMI clustering
echo "Starting UMI clustering..."
for file in $(ls ${out_path}${sample_name}*data2*LTR*.txt); 
do python3 UMI_clustering_hamming_ref.py $file ${file}_hamming_ref_mms${mms}.txt --mismatch_threshold ${mms}; 
done
echo "UMI clustering completed."

#Regroup all the chromosome files together
echo "Merging LTR5 files..."
awk 'FNR==1 && NR!=1 { next } { print }' ${out_path}${sample_name}*data2*_LTR5.txt_hamming_ref_mms${mms}.txt > "${out_path}${sample_name}_merged_LTR5_mms${mms}.txt"
echo "Merging LTR3 files..."
awk 'FNR==1 && NR!=1 { next } { print }' ${out_path}${sample_name}*data2*_LTR3.txt_hamming_ref_mms${mms}.txt > "${out_path}${sample_name}_merged_LTR3_mms${mms}.txt"

#Finish the R script: correct UMI and regroup ShS and IS + merge LTR5 and LTR3 results
echo "Starting second R script..."
Rscript /home/mverneret/work/integration/script/run_IS_script_2.R $R_package_path $sample_name $out_path $maxgapIS $nb_read_per_umi $maxgapShS $mms $win_merge $threshold_raw

echo "Job finished successfully."
