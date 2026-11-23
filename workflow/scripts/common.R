# Shared helpers for the pipeline R scripts. Every script sources this file first.

# Use only the packages of the step's conda environment. R also searches the user library
# (R_LIBS_USER, e.g. ~/Library/R/arm64/4.5/library on macOS) even with --vanilla, and packages
# built there for another R installation crash this R when loaded.
.libPaths(character(0))

# Send each step's output and messages (including errors) to its Snakemake log file.
if (length(snakemake@log) > 0) {
  .log_con <- file(snakemake@log[[1]], open = "wt")
  sink(.log_con, type = "output")
  sink(.log_con, type = "message")
}

source(file.path(dirname(snakemake@scriptdir), "scripts", "tsv.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# vsearch may append ";size=N;" to labels in .uc files
strip_size <- function(x) sub(";size=[0-9]+;?$", "", x)

# Read the S (new cluster) and H (hit) records of a vsearch .uc file.
# Returns query (the sequence) and cluster (its centroid / seed).
read_uc_members <- function(path) {
  read_tsv(path, col_names = FALSE, col_types = cols_only(X1 = "c", X9 = "c", X10 = "c"),
           progress = FALSE) %>%
    filter(X1 %in% c("S", "H")) %>%
    transmute(query = strip_size(X9),
              cluster = if_else(X1 == "S", strip_size(X9), strip_size(X10)))
}

# Sample name of a read label "<sample>-<read number>"
sample_of_read <- function(x) sub("-[0-9]+$", "", x)

write_fasta_simple <- function(names, seqs, path) {
  writeLines(as.vector(rbind(paste0(">", names), seqs)), path)
}
