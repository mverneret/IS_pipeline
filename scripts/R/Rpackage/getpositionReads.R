getpositionReads <- function(file = NULL, randomTag = NULL, LTR = NULL,mapq.val=20){
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(readr))
  
  LTR_file <- merge(file, randomTag, by = "readID", all.x = TRUE)%>%
    mutate(strand=if_else(strand.genome == strand.target, "+","-"))
  
  LTR_file_filteredUMI <- LTR_file[!is.na(LTR_file$Sequence), ]
  print(paste0("Number of reads with UMI : ", nrow(LTR_file_filteredUMI)))
  print(paste0("% of reads with UMI: ", round((nrow(LTR_file_filteredUMI)/ nrow(LTR_file)) * 100,2)))
  
  LTR_file_filtered <- LTR_file_filteredUMI[LTR_file_filteredUMI$mapq.genome>=mapq.val,]
  
  print(paste0("Number of reads conserved with MAPQ >=", mapq.val,": ",nrow(LTR_file_filtered)))
  print(paste0("% of reads conserved with MAPQ >= ", mapq.val, ": ", round((nrow(LTR_file_filtered)/ nrow(LTR_file_filteredUMI)) * 100,2)))
  
  chromosome_list <- split(LTR_file_filtered, LTR_file_filtered$seqnames.genome)
  
  return(chromosome_list)
  
}
