#!/usr/bin/env Rscript
# One-off builder: collect the expert subtype assignments of the legacy job folders into
# resources/subtype_library.tsv (one row per sequence per source).
# Legacy folder conventions are assumed here and nowhere in the pipeline itself.
#
# Every folder of LEGACY_ROOT with exactly one `*subtype_assignments*.xlsx` and one
# `*_final_sequences_to_subtype.fasta` is a source. Sources are written to the library under neutral codes
# (set_01, set_02, ...), so the library does not identify where the data came from. The code to folder name
# map is written to SOURCE_MAP_OUT: keep that file private, outside the repository.
#
# Usage: build_subtype_library.R LEGACY_ROOT PIPELINE_ROOT SOURCE_MAP_OUT [PREFERRED_ASSIGNMENT_REGEX]
#   PREFERRED_ASSIGNMENT_REGEX picks one workbook (matched on the file path) in a folder that has several.

suppressPackageStartupMessages({library(openxlsx); library(digest); library(stringr)})
a <- commandArgs(trailingOnly = TRUE)
if (!length(a) %in% 3:4) stop("usage: build_subtype_library.R LEGACY_ROOT PIPELINE_ROOT SOURCE_MAP_OUT [PREFERRED_ASSIGNMENT_REGEX]")
preferred <- if (length(a) == 4) a[4] else ""
source(file.path(a[2], "workflow", "scripts", "tsv.R"))
source(file.path(a[2], "workflow", "scripts", "subtype_lib.R"))
refs <- readAAStringSet(file.path(a[2], "resources", "ENV_reference_sequence_translations_A-N.fasta"))

jobs <- sort(list.dirs(a[1], recursive = FALSE))
rows <- list()
source_map <- tibble(code = character(), original_source = character())
for (job in jobs) {
  af <- list.files(job, pattern = "subtype_assignments.*\\.xlsx$", full.names = TRUE)
  af <- af[!grepl("~\\$", af)]
  if (nzchar(preferred) && any(grepl(preferred, af))) af <- af[grepl(preferred, af)]
  ff <- list.files(job, pattern = "_final_sequences_to_subtype\\.fasta$", full.names = TRUE)
  if (length(af) != 1 || length(ff) != 1) { message("skip ", basename(job)); next }
  asg <- read.xlsx(af, sheet = 1)
  fa <- readDNAStringSet(ff); names(fa) <- sub(" .*", "", names(fa))
  keep <- asg$new_rep_97_ID %in% names(fa) & !is.na(asg$Subtype) & trimws(asg$Subtype) != ""
  message(basename(job), " (", sprintf("set_%02d", nrow(source_map) + 1), "): ", sum(keep), " of ", nrow(asg), " assignments with sequence")
  asg <- asg[keep, ]
  code <- sprintf("set_%02d", nrow(source_map) + 1)
  source_map <- bind_rows(source_map, tibble(code = code, original_source = basename(job)))
  tr <- translate_best_frame(fa[asg$new_rep_97_ID], refs, cores = 8)
  nt <- toupper(tr$nt_oriented)
  legacy_aa <- if ("Sequence" %in% names(asg)) asg$Sequence else NA_character_
  rows[[code]] <- tibble(
    seq_md5 = vapply(nt, digest, "", algo = "md5", serialize = FALSE), nt_sequence = nt,
    aa_sequence = tr$protein, Subtype = trimws(asg$Subtype),
    Endogenous = if ("Endogenous" %in% names(asg)) as.character(asg$Endogenous) else NA_character_,
    Unsure = as.character(asg$Unsure), Note = as.character(asg$Note), source = code,
    date_added = as.character(Sys.Date()),
    legacy_aa_matches = !is.na(legacy_aa) & legacy_aa == sub("\\*$", "", tr$protein))
}
lib <- bind_rows(rows)
# Only the spelling of non-functional is harmonised ("non-functional", "non_functional"); every other label is kept verbatim.
lib$Subtype <- ifelse(grepl("^non[-_ ]?functional$", tolower(lib$Subtype)), "Non-functional", lib$Subtype)
cat(sprintf("translation reproduces the protein of the legacy sheet for %d of %d rows (%.1f%%)\n",
            sum(lib$legacy_aa_matches), nrow(lib), 100 * mean(lib$legacy_aa_matches)))
lib <- lib %>% select(-legacy_aa_matches) %>% distinct(seq_md5, source, Subtype, .keep_all = TRUE)
write_tsv(lib, file.path(a[2], "resources", "subtype_library.tsv"), na = "")
write_tsv(source_map, a[3], na = "")
cat(nrow(lib), "library rows,", n_distinct(lib$seq_md5), "unique sequences\n")
print(table(lib$Subtype))
conf <- lib %>% group_by(seq_md5) %>% summarise(n_labels = n_distinct(Subtype)) %>% filter(n_labels > 1)
cat(nrow(conf), "sequences carry conflicting labels across sources\n")
