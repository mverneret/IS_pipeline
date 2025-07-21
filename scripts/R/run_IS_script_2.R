# Load Packages 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("GenomicRanges")
#BiocManager::install("Rsamtools")
#install.packages("changepoint")
#install.packages("tidyverse")
#BiocManager::install("GenomicAlignments")
#BiocManager::install("stringdist")

args <- commandArgs(trailingOnly = TRUE)

# Help function
print_help <- function() {
  cat("\nUsage: Rscript run_R_script_2.R <R_package_path> <sample_name> <out.path> <maxgapIS> <mapq.val> <nb_read_per_umi> <maxgapShS> <mms> <win_merge> <threshold.raw>\n")
  cat("\nArguments:\n")
  cat("  R_package_path     Path to the directory containing required R functions\n")
  cat("  sample_name        Name of the sample being processed (ex:'barcode01')\n")
  cat("  out.path           Directory to save output files\n")
  cat("  maxgapIS          Maximum gap allowed between integration sites (IS)\n")
  cat("  nb_read_per_umi   Minimum number of reads per UMI (Unique Molecular Identifier)\n")
  cat("  maxgapShS         Maximum gap allowed between ShearSite (ShS)\n")
  cat("  mms               Mismatches allowed between UMI\n")
  cat("  win_merge         Window size for merging IS\n")
  cat("  threshold.raw     Minimum number of raw reads for filtering IS\n")
  quit(status = 0)
}

# Check if help is requested
if (length(args) == 1 && (args[1] == "--help" || args[1] == "-h")) {
  print_help()
}

# Ensure correct number of arguments
if (length(args) < 9) {
  cat("\nError: Missing arguments!\n")
  print_help()
}

# Load R functions
R_package_path = args[1]
source(paste(R_package_path,"getISposition.R", sep=""))
source(paste(R_package_path,"getIScounts.R", sep=""))
source(paste(R_package_path,"groupGenomicPositions.R", sep=""))
source(paste(R_package_path,"mergeLTRs.R", sep=""))

# Load the variables
sample_name = args[2]
out.path = args[3]

maxgapIS = as.numeric(args[4])
nb_read_per_umi = as.numeric(args[5])
maxgapShS = as.numeric(args[6])
mms = as.numeric(args[7])

win = as.numeric(args[8])
threshold.raw = as.numeric(args[9])

print("1. ------ LTR5 ------")
# import the position files with UMI groups for each chromosome
# Get the IS position of filtered reads grouping reads closer than maxgap after correcting for ShearSite positions and mismatchs in UMIs
LTR="LTR5"
combined_data.LTR5 <- read.table(paste(out.path,sample_name,"_merged_LTR5_mms",mms,".txt", sep=""), header = TRUE, sep = "\t")
print("Get LTR5 IS positions")
LTR5.positions <- getISposition(file=combined_data.LTR5, nb_read_per_umi=nb_read_per_umi, maxgapIS = maxgapIS, maxgapShS=maxgapShS,LTR=LTR)
print("Get LTR5 IS counts")
IS.counted.LTR5 <- getIScounts(LTR5.positions)

print("2. ------ LTR3 ------")
# import the position files with UMI groups for each chromosome
# Get the IS position of filtered reads grouping reads closer than maxgap after correcting for ShearSite positions and mismatchs in UMIs
LTR="LTR3"
combined_data.LTR3 <- read.table(paste(out.path,sample_name,"_merged_LTR3_mms",mms,".txt", sep=""), header = TRUE, sep = "\t")
print("Get LTR3 IS positions")
LTR3.positions <- getISposition(file=combined_data.LTR3, nb_read_per_umi=nb_read_per_umi, maxgapIS = maxgapIS, maxgapShS=maxgapShS, LTR=LTR)
print("Get LTR3 IS counts")
IS.counted.LTR3 <- getIScounts(LTR3.positions)

print("3. ------ Group LTR5-LTR3 results ------")
LTR.positions <- bind_rows(IS.counted.LTR3, IS.counted.LTR5)
LTR.merged <- mergeLTRs(IS=LTR.positions,win=win,threshold.raw = threshold.raw)

print("4. ------ Save Results ------")
# 1. Before count
suppressWarnings(write.table(LTR5.positions, paste(out.path,sample_name, "_positionsreadsLTR5_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

suppressWarnings(write.table(LTR3.positions, paste(out.path,sample_name, "_positionsreadsLTR3_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

# 2. Before merging IS
suppressWarnings(write.table(IS.counted.LTR5, paste(out.path,sample_name, "_countedreadsLTR5_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

suppressWarnings(write.table(IS.counted.LTR3, paste(out.path,sample_name, "_countedreadsLTR3_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

# 3. Before merging LTR5 and LTR3
suppressWarnings(write.table(LTR.positions, paste(out.path,sample_name, "_countedreads_LTR5LTR3_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

# 4. After merging
suppressWarnings(write.table(LTR.merged, paste(out.path,sample_name,"_clonalityResults_mms",mms,"_ShS",maxgapShS,".txt", sep=""),
                             sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))

