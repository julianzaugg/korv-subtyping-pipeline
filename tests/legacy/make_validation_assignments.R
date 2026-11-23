#!/usr/bin/env Rscript
# Validation helper: map a legacy job's expert assignments onto the new run's final
# representatives by sequence content, so the report can be generated without new manual work.
#
# Usage: make_validation_assignments.R LEGACY_DIR LEGACY_PREFIX NEW_DIR NEW_NAME OUT_XLSX

source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "validation_helpers.R"))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5) stop("usage: make_validation_assignments.R LEGACY_DIR LEGACY_PREFIX NEW_DIR NEW_NAME OUT_XLSX")
out <- legacy_truth(a[1], a[2], file.path(a[3], a[4], paste0(a[4], "_final_sequences_to_subtype.tsv")))
# Representatives with no legacy sequence within 5% cannot be given a label
out$Subtype[is.na(out$Subtype)] <- "unassigned_in_legacy"
out$Note <- ifelse(is.na(out$match) | out$match == "exact", NA, paste("validation:", out$match))
cat(sprintf("%d new representatives: %d exact legacy matches, %d nearest, %d unassigned\n", nrow(out),
            sum(out$match == "exact", na.rm = TRUE), sum(grepl("nearest", out$match)), sum(out$Subtype == "unassigned_in_legacy")))
write.xlsx(out %>% select(-match), a[5])
