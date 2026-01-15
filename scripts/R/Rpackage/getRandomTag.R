getRandomTag <- function(fasta_file) {
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(dplyr))
  
  # Define a function to compute the reverse complement of a DNA sequence
  reverse_complement <- function(sequence) {
    # Replace bases with their complements and reverse the sequence
    complement <- chartr("ACGTacgt", "TGCAtgca", sequence)
    paste(rev(strsplit(complement, "")[[1]]), collapse = "")
  }
  
  # Define a function to parse each header line and extract information
  parse_fasta_header <- function(header) {
    # Extract the sequence name, umi_fwd_seq, and umi_rev_seq using regex
    readID <- str_match(header, ">([^;]+)")[, 2]
    umi_fwd_seq <- str_match(header, "umi_fwd_seq=([A-Za-z]*)")[, 2]
    umi_rev_seq <- str_match(header, "umi_rev_seq=([A-Za-z]*)")[, 2]
    
    umi_fwd_seq <- ifelse(umi_fwd_seq %in% c("", "None"), NA, umi_fwd_seq)
    umi_rev_seq <- ifelse(umi_rev_seq %in% c("", "None"), NA, umi_rev_seq)

    # Decide which UMI to use
    if (!is.na(umi_fwd_seq)) {
      umi <- umi_fwd_seq
    } else if (!is.na(umi_rev_seq)) {
      umi <- reverse_complement(umi_rev_seq)
    } else {
      umi <- NA
    }
    
    # Return a data frame with extracted information
    data.frame(
      readID = readID,
      Sequence = umi,
      stringsAsFactors = FALSE
    )
  }
  
  # Read the FASTA file
  fasta_lines <- readLines(fasta_file)
  
  # Filter only header lines
  headers <- fasta_lines[startsWith(fasta_lines, ">")]
  
  # Parse each header and combine into a single data frame
  randomTag <- do.call(rbind, lapply(headers, parse_fasta_header))
  
  # Apply additional sequence modifications
  randomTag$Sequence <- gsub("^TTT|TTT$", "", randomTag$Sequence)
  randomTag$Sequence <- gsub("TT", "", randomTag$Sequence)
  
  return(randomTag)
}
