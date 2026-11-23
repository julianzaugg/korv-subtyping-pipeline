#!/usr/bin/env Rscript
# How good are the Stage 3 suggestions? Compares the suggested subtype of each sequence in a run with
# the expert label of the matching legacy sequence. For an honest test the run must have used a library
# from which the job's own labels were removed (see tests/legacy/make_loo_library.R).
#
# Usage: evaluate_suggestions.R LEGACY_DIR LEGACY_PREFIX NEW_DIR NEW_NAME

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "validation_helpers.R"))
a <- commandArgs(trailingOnly = TRUE)
truth <- legacy_truth(a[1], a[2], file.path(a[3], a[4], paste0(a[4], "_final_sequences_to_subtype.tsv")))
sug <- read.xlsx(file.path(a[3], a[4], paste0(a[4], "_subtype_suggestions.xlsx")), sheet = "suggestions")
d <- sug %>% select(new_rep_97_ID, Total_number_of_reads, suggested_subtype, confidence, basis, flags) %>%
  left_join(truth %>% select(new_rep_97_ID, truth = Subtype, match), by = "new_rep_97_ID") %>%
  filter(!is.na(truth)) %>%
  mutate(sug = ifelse(grepl("^possible new", suggested_subtype), "possible new", suggested_subtype), correct = sug == truth)
cat(sprintf("%d of %d sequences have a legacy label (%d exact sequence matches)\n", nrow(d), nrow(sug), sum(d$match == "exact")))
cat("\nBy confidence tier (sequences / correct / reads covered / reads correct):\n")
print(d %>% group_by(confidence) %>% summarise(n = n(), n_correct = sum(correct), accuracy = round(mean(correct), 3),
                                               reads = sum(Total_number_of_reads), reads_correct = sum(Total_number_of_reads[correct])) %>%
        as.data.frame())
cat(sprintf("\nOverall: %d / %d correct (%.1f%%); read-weighted %.1f%%\n", sum(d$correct), nrow(d), 100 * mean(d$correct),
            100 * sum(d$Total_number_of_reads[d$correct]) / sum(d$Total_number_of_reads)))
cat("\nWrong suggestions (largest first):\n")
print(d %>% filter(!correct) %>% arrange(desc(Total_number_of_reads)) %>% transmute(rep = new_rep_97_ID, reads = Total_number_of_reads, truth, suggested = sug, confidence, flags = substr(flags, 1, 40)) %>% head(25) %>% as.data.frame(), row.names = FALSE)
cat("\nConfusion (rows = expert label, columns = suggestion):\n"); print(table(truth = d$truth, suggested = d$sug))
