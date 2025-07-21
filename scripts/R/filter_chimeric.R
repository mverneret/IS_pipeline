filter_chimeric <- function(minimap2PAF = NULL, targetName = NULL){
  
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(GenomicRanges))
  
  print(paste0("Initial number of reads : ", length(unique(minimap2PAF$qName))))
  
  if(!is_empty(minimap2PAF)){
    
    # 1. Keep primary alignment (tp:A:P)
    print("1. Remove multi-mapped reads with secundary alignments")

    #keep all the primary reads
    minimap2PAF.primaryall <- minimap2PAF %>%
      filter(numberMismatchGaps == "tp:A:P") %>%
      select(-numberMismatchGaps)
    
    #extract the readID of the multi-mapped reads
    minimap2PAF.secundary <- minimap2PAF %>%
      filter(numberMismatchGaps == "tp:A:S") %>%
      select(-numberMismatchGaps)
    
    #filter out the multi-mapped reads
    minimap2PAF.primary <- minimap2PAF %>%
      filter(!(qName %in% unique(minimap2PAF.secundary$qName)))
    
    print(paste0("Number of primary reads: ", length(unique(minimap2PAF.primary$qName))))
    
    # 2. Filter TARGET-HOST reads
    
    if(!is.null(targetName)){
      
      print("2. Remove unwanted target-target and host-host sequences")
      
      paf_filtered <- minimap2PAF.primary %>%
        # 2.1. Add a "lig" field. 1 = TARGET and 0 = HOST.
        mutate(lig = if_else(tName == targetName, 1, 0)) %>%
        group_by(qName) %>%
        # 2.2. Exclude non chimeric reads.
        filter( any(grepl(targetName, tName) ) & !all(tName == targetName)) %>%
        # 2.3. Exclude chimeric reads mapped to several HOST chromosomes.
        filter(length(unique(tName)) < 3) %>%
        arrange(qStart) %>%
        # 2.4. Collapse "lig" information at the read level (10, 01)
        mutate(ligation = paste0(lig , collapse = "")) %>%
        # 2.5. Alignment gaps between substrings
        mutate(mapGap = lead(qStart, default = unique(qLength)) - qEnd) %>%
        ungroup() %>%
        mutate(target = tName %in% targetName) %>%
        select(-lig)
      
      print(paste0("Reads displaying adequat target-genome structure: ", length(unique(paf_filtered$qName))))
      
      print(paste0("Ligation structure Top 10:"))
      print(
        paf_filtered %>% select(qName, ligation) %>%
          distinct() %>%
          group_by(ligation) %>%
          summarise(count = n()) %>%
          arrange(desc(count)) %>%
          dplyr::rename("Ligation type (0 = chrom, 1 = target)" = ligation, "Count" = count) %>%
          head(10)
      )

      
      print("3. Select reads with a '01'|'10' structure")
      
      # 3. Exclude unwanted read structure
      # Add extra formating for GRanges transformation.
      paf_filtered <- paf_filtered %>%
        mutate(width = tEnd - tStart) %>%
        select(tName, tStart, tEnd, qStart, qEnd, width, orientationToRef, ligation, mapGap, MAPQ, qName, qLength, target) %>%
        arrange(qName) %>%
        dplyr::rename("seqnames" = tName,
                      "start" = tStart,
                      "end" = tEnd,
                      "readStart_position" = qStart,
                      "readEnd_position" = qEnd,
                      "strand" = orientationToRef,
                      "readID" = qName) %>%
        dplyr::filter(ligation %in% c("01","10")) %>%
        # Compute the mean and sd gap between substrings.
        group_by(readID) %>%
        mutate(meanGap = mean(abs(mapGap)),
               sdGap = sd(abs(mapGap))) %>%
        ungroup()
      
    } else {
      
      stop('No Target specified! Please provide one')
      
    }
    
  } else{
    
    warning('Empty or Absent minimap2PAF argument!')
    
    emptyFile <- tibble(seqnames = character(),
                        start = numeric(),
                        end = numeric(),
                        readStart_position = numeric(),
                        readEnd_position = numeric(),
                        width = numeric(),
                        strand = character(),
                        ligation = character(),
                        mapGap = character(),
                        MAPQ = character(),
                        readID = character(),
                        qLength = character(),
                        target = character())
    
    paf_filtered <- emptyFile
    
  }
  
  print(paste0("Reads Post-Filtering: ", length(unique(paf_filtered$readID))))
  print(paste0("% Retained : ", 100*round(length(unique(paf_filtered$readID))/length(unique(minimap2PAF$qName)),3)))
  
  return(paf_filtered)
  
}
