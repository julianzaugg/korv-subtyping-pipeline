#!/usr/bin/env Rscript
# Leave-one-job-out library: drop every row whose source matches a pattern, so the suggestions for that
# job cannot see their own expert labels.   Usage: make_loo_library.R LIBRARY.tsv SOURCE_REGEX OUT.tsv
tsv_path <- function() {
  f <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  file.path(dirname(normalizePath(f)), "..", "..", "workflow", "scripts", "tsv.R")
}
source(tsv_path())
suppressPackageStartupMessages(library(dplyr))
a <- commandArgs(trailingOnly = TRUE)
lib <- read_tsv(a[1], col_types = cols(.default = "c"), progress = FALSE)
keep <- lib %>% filter(!grepl(a[2], source))
cat(nrow(lib) - nrow(keep), "rows removed,", nrow(keep), "kept\n")
write_tsv(keep, a[3], na = "")
