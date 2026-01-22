getRandomTag <- function(fasta_file) {
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(dplyr))
  
  # Define a function to compute the reverse complement of a DNA sequence
  reverse_complement <- function(sequence) {
    # Replace bases with their complements and reverse the sequence
    complement <- chartr("ACGT", "TGCA", sequence)
    return(paste(rev(unlist(strsplit(complement, ""))), collapse = ""))
  }
  
  # Define a function to parse each header line and extract information
  parse_fasta_header <- function(header, apply_reverse_complement) {
    # Extract the sequence name, umi_fwd_seq, and umi_rev_seq using regex
    name <- str_match(header, ">([^;]+)")[, 2]
    umi_fwd_seq <- str_match(header, "umi_fwd_seq=([A-Za-z]*)")[, 2]
    umi_rev_seq <- str_match(header, "umi_rev_seq=([A-Za-z]*)")[, 2]
    
    # If apply_reverse_complement is TRUE, take the reverse complement of umi_rev_seq
    umi <- if (apply_reverse_complement) {
      reverse_complement(umi_rev_seq)
    } else {
      umi_fwd_seq
    }
    
    # Return a data frame with extracted information
    data.frame(
      readID = name,
      Sequence = ifelse(umi == "", NA, umi),
      stringsAsFactors = FALSE
    )
  }
  
  # Determine if reverse complement should be applied based on the filename
  apply_reverse_complement <- grepl("LTR3", fasta_file)
  
  # Read the FASTA file
  fasta_lines <- readLines(fasta_file)
  
  # Filter only header lines
  headers <- fasta_lines[startsWith(fasta_lines, ">")]
  
  # Parse each header and combine into a single data frame
  randomTag <- do.call(rbind, lapply(headers, parse_fasta_header, apply_reverse_complement))
  
  # Apply additional sequence modifications
  randomTag$Sequence <- gsub("^TTT|TTT$", "", randomTag$Sequence)
  randomTag$Sequence <- gsub("TT", "", randomTag$Sequence)
  
  return(randomTag)
}

