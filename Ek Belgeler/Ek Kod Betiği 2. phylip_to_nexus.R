#!/usr/bin/env Rscript

# Simple PHYLIP to NEXUS Converter
# This script provides easy-to-use functions for converting PHYLIP to NEXUS

# =============================================================================
# METHOD 1: Using APE package (Recommended - most reliable)
# =============================================================================

convert_phylip_to_nexus_ape <- function(phylip_file, nexus_file = NULL) {
  # Install ape if needed
  if (!require("ape", quietly = TRUE)) {
    install.packages("ape")
    library(ape)
  }
  
  # Set output filename if not provided
  if (is.null(nexus_file)) {
    nexus_file <- sub("\\.phy$|\\.phylip$", ".nex", phylip_file)
  }
  
  # Read PHYLIP and write NEXUS
  alignment <- read.dna(phylip_file, format = "sequential")
  write.nexus.data(alignment, file = nexus_file, format = "dna")
  
  cat("??? Converted:", phylip_file, "???", nexus_file, "\n")
  return(nexus_file)
}

# =============================================================================
# METHOD 2: Manual conversion (no packages required)
# =============================================================================

convert_phylip_to_nexus_manual <- function(phylip_file, nexus_file = NULL) {
  
  if (is.null(nexus_file)) {
    nexus_file <- sub("\\.phy$|\\.phylip$", ".nex", phylip_file)
  }
  
  # Read file
  lines <- readLines(phylip_file)
  lines <- lines[nchar(trimws(lines)) > 0]  # Remove empty lines
  
  # Parse header
  header <- as.integer(strsplit(trimws(lines[1]), "\\s+")[[1]])
  ntax <- header[1]
  nchar <- header[2]
  
  # Parse sequences
  taxa <- character(ntax)
  seqs <- character(ntax)
  
  for (i in 1:ntax) {
    line <- lines[i + 1]
    taxa[i] <- trimws(substr(line, 1, 10))
    seqs[i] <- gsub("\\s+", "", substr(line, 11, nchar(line)))
  }
  
  # Write NEXUS
  cat("#NEXUS\n\n", file = nexus_file)
  cat("BEGIN TAXA;\n", file = nexus_file, append = TRUE)
  cat(paste0("\tDIMENSIONS NTAX=", ntax, ";\n"), file = nexus_file, append = TRUE)
  cat("\tTAXLABELS\n", file = nexus_file, append = TRUE)
  cat(paste0("\t\t", taxa, "\n"), file = nexus_file, append = TRUE)
  cat("\t;\nEND;\n\n", file = nexus_file, append = TRUE)
  
  cat("BEGIN CHARACTERS;\n", file = nexus_file, append = TRUE)
  cat(paste0("\tDIMENSIONS NCHAR=", nchar, ";\n"), file = nexus_file, append = TRUE)
  cat("\tFORMAT MISSING=? GAP=- MATCHCHAR=. DATATYPE=nucleotide;\n", file = nexus_file, append = TRUE)
  cat("\tMATRIX\n", file = nexus_file, append = TRUE)
  
  for (i in 1:ntax) {
    cat(sprintf("%-12s %s\n", taxa[i], seqs[i]), file = nexus_file, append = TRUE)
  }
  
  cat("\t;\nEND;\n", file = nexus_file, append = TRUE)
  
  cat("??? Converted:", phylip_file, "???", nexus_file, "\n")
  return(nexus_file)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

# Example 1: Convert using ape (recommended)
# convert_phylip_to_nexus_ape("alignment.phy")

# Example 2: Convert with custom output name
# convert_phylip_to_nexus_ape("alignment.phy", "output.nex")

# Example 3: Convert without ape package
# convert_phylip_to_nexus_manual("alignment.phy")

# Example 4: Batch convert all PHYLIP files in a directory
batch_convert <- function(directory = ".", method = "ape") {
  phylip_files <- list.files(directory, pattern = "\\.(phy|phylip)$", full.names = TRUE)
  
  if (length(phylip_files) == 0) {
    cat("No PHYLIP files found in", directory, "\n")
    return(invisible(NULL))
  }
  
  cat("Found", length(phylip_files), "PHYLIP file(s)\n")
  
  for (file in phylip_files) {
    if (method == "ape") {
      convert_phylip_to_nexus_ape(file)
    } else {
      convert_phylip_to_nexus_manual(file)
    }
  }
  
  cat("\nBatch conversion complete!\n")
}

# =============================================================================
# COMMAND LINE USAGE
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("\n=== PHYLIP to NEXUS Converter ===\n\n")
    cat("Usage:\n")
    cat("  Rscript convert.R <input.phy> [output.nex]\n\n")
    cat("Examples:\n")
    cat("  Rscript convert.R alignment.phy\n")
    cat("  Rscript convert.R alignment.phy output.nex\n")
    cat("  Rscript convert.R --batch  # Convert all .phy files\n\n")
    quit(status = 1)
  }
  
  if (args[1] == "--batch") {
    batch_convert()
  } else {
    input_file <- args[1]
    output_file <- if (length(args) >= 2) args[2] else NULL
    
    if (!file.exists(input_file)) {
      cat("Error: File not found:", input_file, "\n")
      quit(status = 1)
    }
    
    # Try with ape first, fall back to manual if ape not available
    tryCatch({
      convert_phylip_to_nexus_ape(input_file, output_file)
    }, error = function(e) {
      cat("Note: Using manual conversion (ape package not available)\n")
      convert_phylip_to_nexus_manual(input_file, output_file)
    })
  }
}
#convert_phylip_to_nexus_ape("alignment_manual.phy", nexus_file = "manualnexusalignment.nex")
convert_phylip_to_nexus_manual("alignment_manual.phy", nexus_file = "apenexusalignment.nex")