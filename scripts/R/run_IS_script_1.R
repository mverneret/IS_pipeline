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
  cat("\nUsage: Rscript run_IS_script_1.R <R_package_path> <sample_name> <out.path> <paf.path> <umi.path> <LTR5_targetName> <LTR5_lengthTarget> <LTR3_targetName> <LTR3_lengthTarget> <mapq.val> <assembly>\n")
  cat("\nArguments:\n")
  cat("  R_package_path       Path to the directory containing required R functions\n")
  cat("  sample_name          Name of the sample being processed (e.g., 'barcode01')\n")
  cat("  out.path             Directory to save output files\n")
  cat("  paf.path             Path to the PAF (Pairwise Alignment Format) files\n")
  cat("  umi.path             Path to the UMI (Unique Molecular Identifier) files\n")
  cat("  LTR5_targetName      Name of the target sequence for the 5' LTR\n")
  cat("  LTR5_lengthTarget    Length of the 5' LTR target sequence\n")
  cat("  LTR3_targetName      Name of the target sequence for the 3' LTR\n")
  cat("  LTR3_lengthTarget    Length of the 3' LTR target sequence\n")
  cat("  mapq.val             Minimum mapping quality value\n")
  cat("  assembly             Name of the reference assembly used for mapping (e.g., 'ARS' for *_LTR5_mapped_ARS_SUP.paf)\n")
  quit(status = 0)
}

# Check if help is requested
if (length(args) == 1 && (args[1] == "--help" || args[1] == "-h")) {
  print_help()
}

# Ensure correct number of arguments
if (length(args) < 10) {
  cat("\nError: Missing arguments!\n")
  print_help()
}

# Load R functions
R_package_path = args[1]
source(paste(R_package_path,"readPairwiseAlignmentFile.R", sep=""))
source(paste(R_package_path,"filter_chimeric.R", sep=""))
source(paste(R_package_path,"getBreakPoints.R", sep=""))
source(paste(R_package_path,"getRandomTag.R", sep=""))
source(paste(R_package_path,"getpositionReads.R", sep=""))

# Load the variables
sample_name=args[2]
out.path = args[3]
paf.path = args[4]
assembly = args[11]
umi.path = args[5]

LTR5_targetName = args[6]
LTR5_lengthTarget = as.numeric(args[7])

LTR3_targetName = args[8]
LTR3_lengthTarget = as.numeric(args[9])

mapq.val = as.numeric(args[10])

print("0. ------ Load files and remove reads present in both LTR5 and LTR3 files ------")
PAF_LTR5.path = paste(paf.path,sample_name,"_LTR5_mapped_",assembly,"_SUP.paf", sep="")
PAF_LTR5 <- readPairwiseAlignmentFile(alignFile = PAF_LTR5.path)

PAF_LTR3.path = paste(paf.path,sample_name,"_LTR3_mapped_",assembly,"_SUP.paf", sep="")
PAF_LTR3 <- readPairwiseAlignmentFile(alignFile = PAF_LTR3.path)

# Remove reads present in both LTR5 and LTR3
common_qName <- intersect(PAF_LTR5$qName, PAF_LTR3$qName)
cat("Number of common qName values:", length(common_qName), "\n")
# Remove lines with identical qName
PAF_LTR5 <- PAF_LTR5[!PAF_LTR5$qName %in% common_qName, ]
PAF_LTR3 <- PAF_LTR3[!PAF_LTR3$qName %in% common_qName, ]


print("1. ------ LTR5 ------")
LTR="LTR5"
# 1. Filter The Data To keep Only Chimeric Reads:
# Only the best quality reads are retained here (01, 10).
PAF_LTR5.filter <- filter_chimeric(minimap2PAF = PAF_LTR5, targetName = LTR5_targetName)

# 2. Get The Read Target-Host Junctions And Shear Sites.
PAF_LTR5.breakpoints <- getBreakPoints(PAF = PAF_LTR5.filter, lengthTarget = LTR5_lengthTarget, targetName = LTR5_targetName)

# 3.Get the random Tag/UMI of each read
LTR5_UMI.path=paste(umi.path,sample_name,"_LTR5_UMI.fasta", sep="")
randomTag.LTR5 <- getRandomTag(fasta_file=LTR5_UMI.path)

# 4.Prepare files for UMI clustering
chromosome_list.LTR5 <- getpositionReads(file = PAF_LTR5.breakpoints,randomTag = randomTag.LTR5,LTR = LTR, mapq.val = mapq.val)

for (i in 1:length(chromosome_list.LTR5)) {
  data <- chromosome_list.LTR5[[i]]
  
  # Assuming the second column has the chromosome name
  chromosome_name <- data[1, 2]  # Get the value from the second column (for the first row)
  print(paste("Exporting table for barcode:", sample_name, "and chromosome:", chromosome_name))
  
  # Export the data with the chromosome name from the second column in the filename
  suppressWarnings(write.table(data, 
                               paste(out.path,sample_name, paste("_data2_", chromosome_name,"_",LTR, ".txt", sep=""), sep=""),
                              sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))
}

print("2. ------ LTR3 ------")
LTR="LTR3"
# 1. Filter The Data To keep Only Chimeric Reads:
# Only the best quality reads are retained here (01, 10).
PAF_LTR3.filter <- filter_chimeric(minimap2PAF = PAF_LTR3, targetName = LTR3_targetName)

# 2. Get The Read Target-Host Junctions And Shear Sites.
PAF_LTR3.breakpoints <- getBreakPoints(PAF = PAF_LTR3.filter, lengthTarget = LTR3_lengthTarget, targetName = LTR3_targetName)

# 3.Get the random Tag/UMI of each read
LTR3_UMI.path=paste(umi.path,sample_name,"_LTR3_UMI.fasta", sep="")
randomTag.LTR3 <- getRandomTag(fasta_file=LTR3_UMI.path)

# 4.Prepare files for UMI clustering
chromosome_list.LTR3 <- getpositionReads(file = PAF_LTR3.breakpoints,randomTag = randomTag.LTR3,LTR = LTR, mapq.val = mapq.val)

for (i in 1:length(chromosome_list.LTR3)) {
  data <- chromosome_list.LTR3[[i]]
  
  # Assuming the second column has the chromosome name
  chromosome_name <- data[1, 2]  # Get the value from the second column (for the first row)
  print(paste("Exporting table for barcode:", sample_name, "and chromosome:", chromosome_name))
  
  # Export the data with the chromosome name from the second column in the filename
  suppressWarnings(write.table(data, 
                               paste(out.path, sample_name, paste("_data2_", chromosome_name,"_",LTR, ".txt", sep=""), sep=""),
                               sep="\t", quote=F, row.names=F, col.names=T, append = TRUE))
}

