groupGenomicPositions <- function(genomicPositions = NULL, maxgap = 75, col.dist=NULL){
  
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
  
  # 1. Transform the IS list into GRanges object
  positions.gr <- genomicPositions %>%
    dplyr::rename("seqnames" = seqnames.genome, "start" = !!enquo(col.dist)) %>%
    mutate(end = start,
           strand = "*") %>%
    makeGRangesFromDataFrame(seqnames.field = "seqnames",
                             start.field = "start",
                             end.field = "end",
                             keep.extra.columns = T)
  
  # 2. Group close overlaping IS taken into account a up/down window
  windows <- GenomicRanges::reduce(positions.gr + maxgap)
  
  # 3. Find overlaps
  overlap <- findOverlaps(positions.gr, windows, maxgap = maxgap)
  
  # 4. Reports for each IS the window into which it is falling
  return(subjectHits(overlap))
  
}
