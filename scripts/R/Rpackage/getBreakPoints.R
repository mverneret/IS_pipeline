getBreakPoints <- function(PAF = NULL, targetName = NULL, gapAlignment = NA, distanceToLTR = NA, returnFILTEREDout = FALSE){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))
  
  if(is_empty(PAF)){ stop('Empty or Absent PAF argument!') }
  if(is_empty(targetName)){ stop('Empty or Absent targetName argument! Please provide the TARGET chromosome name') }
  if(isFALSE(all(c("ligation", "mapGap") %in% colnames(PAF)))){ stop('PAF needs to be prepared with filter_chimeric() prior to getBreakPoints()!') }
  
  # 0. OPTION: REMOVE READS WITH AN ALIGNMENT GAP >= gapAlignment
  if(!is.na(gapAlignment)){PAF <- filter(PAF, mapGap < gapAlignment)}
  
  ###################
  #### 1. TARGET ####
  ###################
  if(nrow(PAF) > 100000){ print('This may take some time ...') }
  
  print("1. Summarise TARGET substring information")
  
  targets <- PAF %>%
    select(readID, seqnames, start, end, strand, ligation, readStart_position, readEnd_position, mapGap) %>%
    # 1.1. Filter out GENOME substrings
    filter(seqnames == targetName)%>%
    group_by(readID) %>%
    mutate(strand.target = strand)%>%
    mutate(breakPtsPosition = case_when(
      ligation == "01" ~ "right",
      ligation == "10" ~ "left")) %>%
    ungroup()
  
  ###################
  #### 2. GENOME ####
  ###################
  
  print('2. Summarise GENOME substring information')
  
  genome_concat <- PAF %>%
    select(readID, seqnames, start, end, strand, ligation, mapGap, MAPQ) %>%
    # 2.1. Filter out TARGET substrings
    filter(seqnames != targetName) %>%
    dplyr::rename("seqnames.genome" = seqnames,
                  "start.genome" = start,
                  "end.genome" = end,
                  "strand.genome" = strand,
                  "mapGap.genome"=mapGap,
                  "mapq.genome"=MAPQ) %>%
    # 2.2. Merge with TARGET informations
    merge(targets %>%
            group_by(readID) %>%
            select(readID, strand.target, breakPtsPosition, mapGap) %>%
            distinct(),
          all.x=TRUE,
          by = c("readID", "readID")
    ) %>%
    mutate(minmapgap = pmin(mapGap, mapGap.genome, na.rm = TRUE))%>%
    select(-mapGap,-mapGap.genome)
 
   # 3. Extract TARGET-HOST ShearSites and Breakpoints
  print('3. Get the Integration site position')
    
  # Depending on the LTR position in the read and the genome strand, select the right viral-host breakpoints and shear site
  edges <- mutate(genome_concat, category = paste0(strand.genome, ":", breakPtsPosition)) %>% split(.$category)
    
    if(any(names(edges) %in% "-:left")){
      edges.neg.Left <- mutate(edges[["-:left"]], integrationSite = end.genome, shearSite.genome = start.genome)
    } else{
      edges.neg.Left <- data_frame()
    }
    
    if(any(names(edges) %in% "+:left")){
      edges.pos.Left <- mutate(edges[["+:left"]], integrationSite = start.genome, shearSite.genome = end.genome)
    } else{
      edges.pos.Left <- data_frame()
    }
    
    if(any(names(edges) %in% "-:right")){
      edges.neg.right <- mutate(edges[["-:right"]], integrationSite = start.genome, shearSite.genome = end.genome)
    } else{
      edges.neg.right <- data_frame()
    }
    
    if(any(names(edges) %in% "+:right")){
      edges.pos.right <- mutate(edges[["+:right"]], integrationSite = end.genome, shearSite.genome = start.genome)
    } else{
      edges.pos.right <- data_frame()
    }

  read_breakPoints <- bind_rows(edges.neg.Left, edges.pos.Left, edges.neg.right, edges.pos.right) %>%
    select(readID, seqnames.genome, shearSite.genome, integrationSite, strand.target, strand.genome, ligation, minmapgap, mapq.genome)
  
  return(read_breakPoints)

print('5. Done! ')
}
