#!/usr/bin/env Rscript
# Tier A validation of Stage 4: compare the new final workbook with a legacy *_modified workbook.
# Samples are matched on the text before the first "_" (legacy LIMS ID), numbers on exact equality
# for counts and 0.01 tolerance for relative abundances.
#
# Usage: validate_report.R LEGACY_MODIFIED_XLSX NEW_FINAL_XLSX [NEW_SHEET] [LEGACY_SHEET] [A|B]
# Mode A (default): exact counts. Mode B (new random subsample from raw reads): relative abundance
# within 2 percentage points, total KoRV reads within 5%.

suppressPackageStartupMessages({library(dplyr); library(openxlsx)})
a <- commandArgs(trailingOnly = TRUE)
new_sheet <- if (length(a) >= 3) a[3] else "counts_and_abundances"
leg_sheet <- if (length(a) >= 4) a[4] else "counts_and_abundances"
mode <- if (length(a) >= 5) a[5] else "A"

read_ca <- function(path, sheet) {
  m <- read.xlsx(path, sheet = sheet, colNames = FALSE, skipEmptyRows = FALSE)
  h <- which(apply(m, 1, function(r) any(r == "Total_KoRV_reads", na.rm = TRUE)))[1]  # layouts differ between jobs
  if (is.na(h)) stop("no Total_KoRV_reads header found in ", path)
  hdr <- as.character(unlist(m[h, ])); m <- m[-seq_len(h), , drop = FALSE]
  m <- m[!is.na(m[[1]]) & !grepl("^Number of", m[[1]]), , drop = FALSE]
  tot <- which(hdr == "Total_KoRV_reads"); k <- (length(hdr) - tot) / 2
  cats <- hdr[(tot + 1):(tot + k)]
  list(id = sub("_.*", "", as.character(m[[1]])), total = as.numeric(m[[tot]]),
       counts = sapply(seq_len(k), function(i) as.numeric(m[[tot + i]])) %>% matrix(ncol = k, dimnames = list(NULL, cats)),
       rel = sapply(seq_len(k), function(i) as.numeric(m[[tot + k + i]])) %>% matrix(ncol = k, dimnames = list(NULL, cats)))
}
leg <- read_ca(a[1], leg_sheet); new <- read_ca(a[2], new_sheet)
ok <- TRUE
check <- function(label, pass, detail = "") { cat(sprintf("[%s] %s %s\n", if (pass) "PASS" else "FAIL", label, detail)); if (!pass) ok <<- FALSE }

check("same samples", setequal(leg$id, new$id), sprintf("(legacy %d, new %d)", length(leg$id), length(new$id)))
check("same subtype categories", setequal(colnames(leg$counts), colnames(new$counts)),
      sprintf("(legacy: %s | new: %s)", paste(colnames(leg$counts), collapse = ","), paste(colnames(new$counts), collapse = ",")))
common <- intersect(leg$id, new$id); cats <- intersect(colnames(leg$counts), colnames(new$counts))
li <- match(common, leg$id); ni <- match(common, new$id)
dr <- abs(leg$rel[li, cats, drop = FALSE] - new$rel[ni, cats, drop = FALSE])
if (mode == "A") {
  check("Total_KoRV_reads identical", all(leg$total[li] == new$total[ni]),
        sprintf("(%d of %d samples differ)", sum(leg$total[li] != new$total[ni]), length(common)))
  dc <- abs(leg$counts[li, cats, drop = FALSE] - new$counts[ni, cats, drop = FALSE])
  check("read counts per subtype identical", all(dc == 0, na.rm = TRUE), sprintf("(%d cells differ, max diff %g)", sum(dc > 0, na.rm = TRUE), max(dc, na.rm = TRUE)))
  check("relative abundance within 0.01", all(dr <= 0.01, na.rm = TRUE), sprintf("(max diff %g)", max(dr, na.rm = TRUE)))
} else {
  ratio <- new$total[ni] / leg$total[li]
  check("Total_KoRV_reads within 5%", all(abs(ratio - 1) <= 0.05), sprintf("(ratio range %.3f-%.3f)", min(ratio), max(ratio)))
  check("relative abundance within 2 percentage points", all(dr <= 2, na.rm = TRUE), sprintf("(max diff %g)", max(dr, na.rm = TRUE)))
  cmp <- data.frame(sample = leg$id[li], total_legacy = leg$total[li], total_new = new$total[ni],
                    max_abs_diff_pct = apply(dr, 1, max, na.rm = TRUE))
  cmp$max_abs_diff_pct <- round(cmp$max_abs_diff_pct, 2)
  print(cmp, row.names = FALSE)
}
cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
quit(status = if (ok) 0 else 1)
