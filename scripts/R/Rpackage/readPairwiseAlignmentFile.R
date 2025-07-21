readPairwiseAlignmentFile <- function(alignFile = NULL){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))
  
  # 1. LOAD FILE
  print(paste0("Read Pairwise mApping File located at: ", alignFile))
  
  alignMinimap2 <- suppressWarnings(read_delim(alignFile, delim="\t", col_names = c("qName",
                                                                                    "qLength",
                                                                                    "qStart",
                                                                                    "qEnd",
                                                                                    "orientationToRef",
                                                                                    "tName",
                                                                                    "tLength",
                                                                                    "tStart",
                                                                                    "tEnd",
                                                                                    "matchingBases",
                                                                                    "matchingBasesIncludingGaps",
                                                                                    "MAPQ",
                                                                                    "numberMismatchGaps",
                                                                                    "DPscoreMaxScoring",
                                                                                    "DPalignmentScore",
                                                                                    "numberAmbiguousBases",
                                                                                    "alignmentType",
                                                                                    "minimizers",
                                                                                    "X19", "X20","X21", "X22", "X23"),
                                               col_types =
                                                 cols_only(
                                                   qName = col_character(),
                                                   qLength = col_integer(),
                                                   qStart = col_integer(),
                                                   qEnd = col_integer(),
                                                   orientationToRef = col_character(),
                                                   tName = col_character(),
                                                   tLength = col_integer(),
                                                   tStart = col_integer(),
                                                   tEnd = col_integer(),
                                                   matchingBases = col_integer(),
                                                   matchingBasesIncludingGaps = col_integer(),
                                                   MAPQ = col_integer(),
                                                   numberMismatchGaps = col_character(),
                                                   DPscoreMaxScoring = col_character(),
                                                   DPalignmentScore = col_character(),
                                                   numberAmbiguousBases = col_character(),
                                                   alignmentType = col_character(),
                                                   minimizers = col_character()
                                                 )
  )
  )
  
  print(paste0("Number of reads loaded : ", length(unique(alignMinimap2$qName))))
  
  return(alignMinimap2)
  
}
