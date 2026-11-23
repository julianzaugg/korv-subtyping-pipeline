#!/usr/bin/env Rscript
# Tier A validation of Stage 2: compare new clustering outputs with a legacy (QIIME) job.
# Everything is compared by sequence content, never by read ID, because the read chosen
# as centroid can legitimately differ.
#
# Usage: validate_cluster.R LEGACY_DIR LEGACY_PREFIX NEW_DIR NEW_NAME
#   e.g. validate_cluster.R reference_jobs/JOB_A JOB_A results/job_a JOB_A

tsv_path <- function() {
  f <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))
  file.path(dirname(normalizePath(f)), "..", "..", "workflow", "scripts", "tsv.R")
}
source(tsv_path())
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(openxlsx); library(Biostrings)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: validate_cluster.R LEGACY_DIR LEGACY_PREFIX NEW_DIR NEW_NAME")
ldir <- args[1]; lp <- args[2]; ndir <- args[3]; nn <- args[4]
L <- function(suffix) file.path(ldir, paste0(lp, "_", suffix))
N <- function(suffix) file.path(ndir, nn, paste0(nn, "_", suffix))

ok <- TRUE
check <- function(label, pass, detail = "") {
  cat(sprintf("[%s] %s %s\n", if (pass) "PASS" else "FAIL", label, detail))
  if (!pass) ok <<- FALSE
}
fa <- function(path) {
  x <- readDNAStringSet(path)
  tibble(id = sub(";size=[0-9]+;?$", "", names(x)), seq = toupper(as.character(x)))
}

# ---- 100% table, keyed by sequence ---------------------------------------------------
leg_seq <- fa(file.path(L("rep_seq_100.fasta"), "dna-sequences.fasta"))
leg_tab <- read.delim(L("otu_table_100.tsv"), skip = 1, comment.char = "", check.names = FALSE, sep = "\t")
names(leg_tab)[1] <- "id"
leg_long <- leg_tab %>%
  pivot_longer(-id, names_to = "sample", values_to = "count") %>% filter(count > 0) %>%
  left_join(leg_seq, by = "id") %>% group_by(seq, sample) %>% summarise(count = sum(count), .groups = "drop")

new_seq <- read_tsv(file.path(ndir, nn, "work", "otu_100_sequences.tsv"), col_types = "cc", progress = FALSE)
new_long <- read_tsv(file.path(ndir, nn, "work", "otu_table_100_long.tsv.gz"), col_types = "ccd", progress = FALSE) %>%
  filter(count > 0) %>% left_join(new_seq, by = "otu") %>%
  group_by(seq = sequence, sample) %>% summarise(count = sum(count), .groups = "drop")

check("100%: total reads", sum(leg_long$count) == sum(new_long$count),
      sprintf("(legacy %d, new %d)", sum(leg_long$count), sum(new_long$count)))
check("100%: number of OTUs", n_distinct(leg_long$seq) == n_distinct(new_long$seq),
      sprintf("(legacy %d, new %d)", n_distinct(leg_long$seq), n_distinct(new_long$seq)))
check("100%: identical set of representative sequences", setequal(leg_long$seq, new_long$seq))
j <- full_join(leg_long, new_long, by = c("seq", "sample"), suffix = c(".leg", ".new")) %>%
  mutate(across(starts_with("count"), ~ coalesce(.x, 0)))
check("100%: identical count per sequence per sample", all(j$count.leg == j$count.new),
      sprintf("(%d of %d cells differ)", sum(j$count.leg != j$count.new), nrow(j)))

# ---- filtered table ------------------------------------------------------------------
leg_f <- read.csv(L("otu_table_100_filtered_counts.csv"))
new_f <- read.csv(N("otu_table_100_filtered_counts.csv"))
check("filter: per-sample pre/post counts identical",
      isTRUE(all.equal(leg_f[order(leg_f$Sample), ], new_f[order(new_f$Sample), ], check.attributes = FALSE)))
leg_keep <- fa(L("rep_sequences_100_greater_than_2.fasta"))
new_keep <- fa(N("rep_sequences_100_greater_than_2.fasta"))
check("filter: identical set of kept sequences", setequal(leg_keep$seq, new_keep$seq),
      sprintf("(legacy %d, new %d)", n_distinct(leg_keep$seq), n_distinct(new_keep$seq)))

# ---- 97% partition -------------------------------------------------------------------
leg_m <- read.xlsx(L("otu_table_100_with_reps.xlsx"), sheet = "OTU_table_100_filt_with_reps") %>%
  select(rep_100_ID = 1, leg97 = rep_97_ID)
new_m <- read.xlsx(N("otu_table_100_with_reps.xlsx"), sheet = "OTU_table_100_filt_with_reps") %>%
  select(rep_100_ID = 1, new97 = rep_97_ID)
leg_m <- leg_m %>% left_join(leg_keep, by = c("rep_100_ID" = "id"))
new_m <- new_m %>% left_join(new_keep, by = c("rep_100_ID" = "id"))
pm <- inner_join(leg_m %>% select(seq, leg97), new_m %>% select(seq, new97), by = "seq")
n_leg <- n_distinct(pm$leg97); n_new <- n_distinct(pm$new97); n_pair <- nrow(distinct(pm, leg97, new97))
check("97%: number of clusters", n_leg == n_new, sprintf("(legacy %d, new %d)", n_leg, n_new))
check("97%: identical partition of 100% reps", n_leg == n_new && n_pair == n_leg,
      sprintf("(%d distinct legacy/new cluster pairs)", n_pair))

# ---- final sequences -----------------------------------------------------------------
leg_fin <- fa(L("final_sequences_to_subtype.fasta"))
new_fin <- fa(N("final_sequences_to_subtype.fasta"))
check("final: identical set of sequences to subtype", setequal(leg_fin$seq, new_fin$seq),
      sprintf("(legacy %d, new %d; only legacy %d, only new %d)", nrow(leg_fin), nrow(new_fin),
              length(setdiff(leg_fin$seq, new_fin$seq)), length(setdiff(new_fin$seq, leg_fin$seq))))

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
