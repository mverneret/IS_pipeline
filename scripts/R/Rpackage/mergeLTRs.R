mergeLTRs <- function(IS = NULL, win = 25, threshold.raw = 3, threshold.ShS = 2){
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(readr))
  
  # 1. Group each IS located in a +/- 25 bp
  # Select the position of LTR3 for each group
  # raw/filtered: Take the max between the LTRs raw/filtered
  # raw_recall/filtered_recall: Take the max between the LTRs raw_recall/filtered_recall
  IS.final <- IS %>%
    mutate(LTRgroup = groupGenomicPositions(., maxgap = win, col.dist=integrationSite)) %>%
    group_by(LTRgroup,LTR) %>%
    filter(filtered.corShSUMIcor == max(filtered.corShSUMIcor)) %>%
    ungroup() %>%
    group_by(LTRgroup) %>%
    mutate(chr = ifelse(any(LTR %in% "LTR5"), seqnames.genome[LTR == "LTR5"], seqnames.genome),
           exactPosition = ifelse(any(LTR %in% "LTR5") & any(LTR %in% "LTR3"), integrationSite[which.max(filtered.corShSUMIcor)], integrationSite),
           strand = ifelse(any(LTR %in% "LTR5") & any(LTR %in% "LTR3"),
                           ifelse(strand[LTR == "LTR5"] == strand[LTR == "LTR3"],strand[LTR == "LTR5"],"*"),
                           ifelse(any(LTR %in% "LTR5"), strand[LTR == "LTR5"], strand)),
           minPosition.LTR5 = ifelse(any(LTR %in% "LTR5"), minPosition[LTR == "LTR5"], 0),
           maxPosition.LTR5 = ifelse(any(LTR %in% "LTR5"), maxPosition[LTR == "LTR5"], 0),
           minPosition.LTR3 = ifelse(any(LTR %in% "LTR3"), minPosition[LTR == "LTR3"], 0),
           maxPosition.LTR3 = ifelse(any(LTR %in% "LTR3"), maxPosition[LTR == "LTR3"], 0),
           #
           raw.ALL.LTR5 = ifelse(any(LTR %in% "LTR5"), max(raw.ALL[LTR == "LTR5"]), 0),
           raw.ALL.LTR3 = ifelse(any(LTR %in% "LTR3"), max(raw.ALL[LTR == "LTR3"]), 0),
           #
           filtered.rawUMI.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.rawUMI[LTR == "LTR5"]), 0),
           filtered.corUMI.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.corUMI[LTR == "LTR5"]), 0),
           filtered.ShS.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.ShS[LTR == "LTR5"]), 0),
           filtered.corShS.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.corShS[LTR == "LTR5"]), 0),
           filtered.corShSUMIraw.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.corShSUMIraw[LTR == "LTR5"]), 0),
           filtered.corShSUMIcor.LTR5 = ifelse(any(LTR %in% "LTR5"), max(filtered.corShSUMIcor[LTR == "LTR5"]), 0),
           #
           filtered.rawUMI.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.rawUMI[LTR == "LTR3"]), 0),
           filtered.corUMI.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.corUMI[LTR == "LTR3"]), 0),
           filtered.ShS.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.ShS[LTR == "LTR3"]), 0),
           filtered.corShS.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.corShS[LTR == "LTR3"]), 0),
           filtered.corShSUMIraw.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.corShSUMIraw[LTR == "LTR3"]), 0),
           filtered.corShSUMIcor.LTR3 = ifelse(any(LTR %in% "LTR3"), max(filtered.corShSUMIcor[LTR == "LTR3"]), 0),
           #
           raw.ALL.max=max(raw.ALL.LTR5,raw.ALL.LTR3),
           filtered.rawUMI.max = max(filtered.rawUMI.LTR5, filtered.rawUMI.LTR3),
           filtered.corUMI.max = max(filtered.corUMI.LTR5, filtered.corUMI.LTR3),
           filtered.ShS.max = max(filtered.ShS.LTR5, filtered.ShS.LTR3),
           filtered.corShS.max = max(filtered.corShS.LTR5, filtered.corShS.LTR3),
           filtered.corShSUMIraw.max = max(filtered.corShSUMIraw.LTR5, filtered.corShSUMIraw.LTR3),
           filtered.corShSUMIcor.max = max(filtered.corShSUMIcor.LTR5, filtered.corShSUMIcor.LTR3),
           #
           LTR = ifelse(all(LTR == "LTR5"), "LTR5", ifelse(all(LTR == "LTR3"), "LTR3", "LTR5LTR3"))) %>%
    ungroup() %>%
    select(-raw.ALL,
           -filtered.rawUMI,
           -filtered.corUMI,
           -filtered.ShS,
           -filtered.corShS,
           -filtered.corShSUMIraw,
           -filtered.corShSUMIcor,
           -LTRgroup) %>%
    distinct() %>%
    # Final rearrangement:
    select(
      chr,
      strand,
      exactPosition,
      minPosition.LTR5,
      maxPosition.LTR5,
      minPosition.LTR3,
      maxPosition.LTR3,
      LTR,
      raw.ALL.max,
      filtered.corShS.max,
      filtered.corShSUMIcor.max)
  
    #A modifier pour filtrer les site d'intégration et les pourcentages de quantif à afficher 
  IS.final2 <- IS.final %>%
    filter(LTR == "LTR5LTR3" | (raw.ALL.max >= threshold.raw | filtered.corShS.max >= threshold.ShS))%>%
    arrange(desc(filtered.corShSUMIcor.max))%>%
    distinct()
  
  IS.final3 <- IS.final2 %>%
    mutate(perc.corShSUMIcor = round(100*(filtered.corShSUMIcor.max/sum(filtered.corShSUMIcor.max)), 3))

  return(IS.final3)
  
}


