getISposition <- function(file, nb_read_per_umi=2, maxgapIS = 15, maxgapShS=2, LTR){

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(readr))

  #filter the file to keep only UMI (UMI_group so corrected UMI) supported by a minimum number of reads
  final_result <- file %>%
    group_by(UMI_group) %>%
    filter(n() >= nb_read_per_umi) %>%
    ungroup()

  print(paste0("Number of reads conserved with UMI rattached to:",nb_read_per_umi,"reads:",nrow(final_result)))

  IS.groupedIS <- final_result %>%
    mutate(readGroupsIS = groupGenomicPositions(., maxgap = maxgapIS, col.dist=integrationSite)) %>%
    group_by(readGroupsIS)

  IS.groupedShS <- final_result %>%
    mutate(readGroupsShS = groupGenomicPositions(., maxgap = maxgapShS, col.dist=shearSite.genome)) %>%
    group_by(readGroupsShS)

  IS.grouped <- merge(IS.groupedIS,IS.groupedShS) %>%
    mutate(LTR=LTR)

  return(IS.grouped)

}
