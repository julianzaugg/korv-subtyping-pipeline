tsv_path <- function() {
  f <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  file.path(dirname(normalizePath(f)), "..", "..", "workflow", "scripts", "tsv.R")
}
source(tsv_path())
# Shared helpers for the validation scripts.
suppressPackageStartupMessages({library(dplyr); library(openxlsx); library(Biostrings)})

# Legacy expert labels for the sequences of a new run, matched by sequence content.
# Exact nucleotide match first (either strand); otherwise the nearest legacy sequence within
# max_rel_dist edit distance, flagged in `match`. Returns the new_rep_97_ID with Subtype etc.
legacy_truth <- function(legacy_dir, legacy_prefix, new_tsv, max_rel_dist = 0.05) {
  L <- function(s) file.path(legacy_dir, paste0(legacy_prefix, "_", s))
  fa <- function(path) { x <- readDNAStringSet(path); tibble(id = names(x), seq = toupper(as.character(x))) }
  legacy_seq <- fa(L("final_sequences_to_subtype.fasta"))
  assign_files <- Sys.glob(file.path(legacy_dir, "*subtype_assignments*.xlsx")); assign_files <- assign_files[!grepl("~\\$", assign_files)]
  # when a legacy folder holds several workbooks, KORV_PREFERRED_ASSIGNMENT_PATTERN (a regex on the file path) picks one
  preferred <- Sys.getenv("KORV_PREFERRED_ASSIGNMENT_PATTERN", "")
  if (nzchar(preferred) && any(grepl(preferred, assign_files))) assign_files <- assign_files[grepl(preferred, assign_files)]
  legacy <- read.xlsx(assign_files[1], sheet = 1) %>% select(-any_of(c("LIMS_ID", "Sample_ID", "Sequence", "sample_id")))
  legacy <- legacy %>% left_join(legacy_seq, by = c("new_rep_97_ID" = "id")) %>% filter(!is.na(seq))
  legacy$Subtype <- ifelse(grepl("^non[-_ ]?functional$", tolower(trimws(legacy$Subtype))), "Non-functional", trimws(legacy$Subtype))
  ref <- legacy %>% distinct(seq, .keep_all = TRUE)
  rc <- function(x) as.character(reverseComplement(DNAStringSet(x)))
  new <- read_tsv(new_tsv, col_types = cols(.default = "c"), progress = FALSE)
  out <- new %>% select(new_rep_97_ID, seq_md5, sequence) %>%
    left_join(ref %>% select(-new_rep_97_ID), by = c("sequence" = "seq")) %>% mutate(match = ifelse(is.na(Subtype), NA_character_, "exact"))
  miss <- which(is.na(out$Subtype))
  if (length(miss)) {
    d <- adist(out$sequence[miss], ref$seq)
    for (i in seq_along(miss)) {
      j <- which.min(d[i, ]); rel <- d[i, j] / nchar(out$sequence[miss[i]])
      if (rel <= max_rel_dist) {
        for (col in intersect(c("Subtype", "Endogenous"), names(ref))) out[[col]][miss[i]] <- ref[[col]][j]
        out$match[miss[i]] <- sprintf("nearest legacy sequence, %d edits", d[i, j])
      }
    }
  }
  out
}
