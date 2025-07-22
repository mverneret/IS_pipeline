#!/bin/bash
#Lauch on the genobioinfo cluster to use direclty R script

#SBATCH --job-name=R_clonality_random
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH -p workq
#SBATCH --output=run_IS_random.out

module load statistics/R/4.3.1

sample_name="random"
R_package_path="/home/mverneret/work/integration/script/Rpackages/"
out_path="/home/mverneret/work/integration/results_simulation/R_clonality_random/"
input_paf_path="/home/mverneret/work/integration/results_simulation/minimap2/"
input_umi_path="/home/mverneret/work/integration/results_simulation/extract_UMI/"
assembly="ARS12"

###LTR5###
targetName_LTR5="LTR5_sens"
lengthTarget_LTR5=63

###LTR3###
targetName_LTR3="LTR3_sens"
lengthTarget_LTR3=188

maxgapIS=15
mapq_val=20
nb_read_per_umi=1
maxgapShS=10
mms=0
threshold_raw=1
win_merge=1


#Run the first part of the R script : keep chimeric reads, get junctions, add UMI
Rscript /home/mverneret/work/integration/script/run_IS_script_1.R $R_package_path $sample_name $out_path $input_paf_path $input_umi_path $targetName_LTR5 $lengthTarget_LTR5 $targetName_LTR3 $lengthTarget_LTR3 $mapq_val $assembly

#UMI clustering
for file in $(ls ${out_path}${sample_name}*data2*LTR*.txt); 
do python3 /home/mverneret/work/integration/script/UMI_clustering_hamming_ref.py $file ${file}_hamming_ref_mms${mms}.txt --mismatch_threshold ${mms}; 
done
#Regroup all the chromosome files together
awk 'FNR==1 && NR!=1 { next } { print }' ${out_path}${sample_name}*data2*_LTR5.txt_hamming_ref_mms${mms}.txt > "${out_path}${sample_name}_merged_LTR5_mms${mms}.txt"
awk 'FNR==1 && NR!=1 { next } { print }' ${out_path}${sample_name}*data2*_LTR3.txt_hamming_ref_mms${mms}.txt > "${out_path}${sample_name}_merged_LTR3_mms${mms}.txt"

#Finish the R script: correct UMI and regroup ShS and IS + merge LTR5 and LTR3 results
Rscript /home/mverneret/work/integration/script/run_IS_script_2.R $R_package_path $sample_name $out_path $maxgapIS $nb_read_per_umi $maxgapShS $mms $win_merge $threshold_raw

