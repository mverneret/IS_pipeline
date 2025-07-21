getIScounts <- function(file){
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(readr))
  
  IS.counted <- file %>%
    group_by(readGroupsIS) %>%
    mutate(minPosition = min(integrationSite),
           maxPosition= max(integrationSite),
           raw.ALL = n(),
           filtered.rawUMI = length(unique(paste(Sequence))),
           filtered.corUMI = length(unique(paste(UMI_group))),
           filtered.ShS = length(unique(paste(shearSite.genome))),
           filtered.corShS = length(unique(paste(readGroupsShS))),
           filtered.corShSUMIraw = length(unique(paste(readGroupsShS, Sequence))),
           filtered.corShSUMIcor = length(unique(paste(readGroupsShS, UMI_group))),
           filtered.ShSUMIraw = length(unique(paste(shearSite.genome, Sequence))),
           filtered.ShSUMIcor = length(unique(paste(shearSite.genome, UMI_group)))
           
    ) %>%
    dplyr::filter(filtered.corShS == max(filtered.corShS)) %>%
    arrange(desc(filtered.corShS)) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    select(-readID, -ligation, -minmapgap, -mapq.genome, -strand.target, -strand.genome, -Sequence,-UMI_group,-readGroupsIS, -readGroupsShS)
  
  print(paste0("Number of reads without PCR duplicates (mms=0) : ", sum(IS.counted$filtered.rawUMI)))
  print(paste0("Number of reads without PCR duplicates correcting UMIs : ", sum(IS.counted$filtered.corUMI)))
  print(paste0("Number of corrected shearsite : ", sum(IS.counted$filtered.corShS)))
  print(paste0("Number of sister cells (ShS and UMI corrected): ", sum(IS.counted$filtered.corShSUMIcor)))
  print(paste0("Number of integration sites : ", nrow(IS.counted)))

  return(IS.counted)
}
