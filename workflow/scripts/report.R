source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
suppressPackageStartupMessages(library(openxlsx))

# The assignments file is the completed suggestions workbook (.xlsx, also written by the review page)
# or a table with the same columns as .tsv or .csv.
read_assignments <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) as.data.frame(read_tsv(path, col_types = cols(.default = "c"), progress = FALSE), check.names = FALSE)
  else if (ext == "csv") read.csv(path, colClasses = "character", check.names = FALSE, na.strings = "")
  else read.xlsx(path, sheet = 1, check.names = FALSE)
}
assignments <- read_assignments(snakemake@input[["assignments"]])
collapsed <- read.xlsx(snakemake@input[["collapsed"]], sheet = "OTU_table_100_agg_per_rep_97", check.names = FALSE)
final_tsv <- read_tsv(snakemake@input[["final"]], col_types = cols(.default = "c"), progress = FALSE)
filtered_counts <- read.csv(snakemake@input[["filtered_counts"]], check.names = FALSE)
sample_map <- if ("sample_map" %in% names(snakemake@input)) {
  read_tsv(snakemake@input[["sample_map"]], col_types = cols(.default = "c"), progress = FALSE)
} else NULL
stats_file <- if ("read_stats" %in% names(snakemake@input)) snakemake@input[["read_stats"]] else NULL

# ---- assignments: required columns, every representative assigned ---------------------
if (!"Subtype" %in% names(assignments)) stop("assignments file has no 'Subtype' column")
key <- if ("new_rep_97_ID" %in% names(assignments)) "new_rep_97_ID" else
  if ("seq_md5" %in% names(assignments)) "seq_md5" else stop("assignments need new_rep_97_ID or seq_md5")
assignments <- assignments %>% mutate(across(everything(), as.character))
assignments$Subtype <- trimws(assignments$Subtype)
assign_key <- assignments %>% filter(!is.na(Subtype), Subtype != "") %>% distinct(.data[[key]], .keep_all = TRUE)
final_key <- if (key == "seq_md5") "seq_md5" else "new_rep_97_ID"
stale <- setdiff(assignments[[key]][!is.na(assignments[[key]]) & assignments[[key]] != ""], final_tsv[[final_key]])
if (length(stale) > 0) warning(length(stale), " assignment row(s) are not in this run (is the file from another analysis?), e.g. ", paste(head(stale, 3), collapse = ", "))
missing <- final_tsv %>% filter(!(.data[[final_key]] %in% assign_key[[key]]))
if (nrow(missing) > 0) {
  stop(nrow(missing), " representative sequences have no Subtype in the assignments file, e.g. ",
       paste(head(missing$new_rep_97_ID, 5), collapse = ", "), ". Complete the assignments and re-run.")
}
extra_cols <- intersect(snakemake@params[["extra_categories"]], names(assignments))

# ---- long counts per representative and sample -----------------------------------------
sample_cols <- setdiff(names(collapsed), c("new_rep_97_ID", "Cluster_ID", "N_samples_present", "Total_number_of_reads"))
# report the samples in sample sheet order (any not in the sheet, e.g. from a legacy FASTA, at the end)
if (!is.null(sample_map)) sample_cols <- c(intersect(sample_map$internal_id, sample_cols), setdiff(sample_cols, sample_map$internal_id))
display_name <- function(internal) {
  if (is.null(sample_map)) return(internal)
  sample_map$sample_id[match(internal, sample_map$internal_id)]
}
rep_to_category <- final_tsv %>%
  select(new_rep_97_ID, all_of(final_key)) %>%
  left_join(assign_key %>% select(all_of(key), Subtype, any_of(extra_cols)), by = setNames(key, final_key))
if ("Endogenous" %in% extra_cols) {
  # Tidy the free-text Endogenous values: trim spaces, treat spellings that differ only in case as one
  # (the most frequent spelling wins). If the column is entirely blank it is ignored, otherwise blanks
  # become "unassigned" so every representative falls in a category.
  endo <- trimws(as.character(rep_to_category$Endogenous)); endo[is.na(endo)] <- ""
  if (all(endo == "")) {
    message("Endogenous column is empty: only the Subtype summary is written")
    extra_cols <- setdiff(extra_cols, "Endogenous")
  } else {
    lower <- tolower(endo)
    canon <- tapply(endo[endo != ""], lower[endo != ""], function(v) names(sort(table(v), decreasing = TRUE))[1])
    endo <- ifelse(endo == "", "unassigned", unname(canon[lower]))
    if (any(lower == "")) message(sum(lower == ""), " representatives have a blank Endogenous value and are reported as 'unassigned'")
    rep_to_category$Endogenous <- endo
    rep_to_category$Subtype_Endogenous <- paste(rep_to_category$Subtype, endo, sep = "_")
  }
}
long <- collapsed %>%
  select(new_rep_97_ID, all_of(sample_cols)) %>%
  pivot_longer(-new_rep_97_ID, names_to = "internal_id", values_to = "Read_count") %>%
  left_join(rep_to_category, by = "new_rep_97_ID")

# ---- counts and abundances for one category column --------------------------------------
summarise_category <- function(category) {
  per <- long %>%
    group_by(internal_id, category = .data[[category]]) %>%
    summarise(Read_count = sum(Read_count), .groups = "drop")
  cats <- sort(unique(per$category))
  counts <- per %>% pivot_wider(names_from = category, values_from = Read_count, values_fill = 0) %>%
    select(internal_id, all_of(cats))
  counts <- counts[match(sample_cols, counts$internal_id), ]
  total <- rowSums(counts[, cats, drop = FALSE])
  rel <- sweep(as.matrix(counts[, cats, drop = FALSE]), 1, total, "/")
  rel <- round(round(rel, 4) * 100, 2)
  rel[total == 0, ] <- NA
  out <- tibble(sample_id = display_name(counts$internal_id))
  if (!is.null(sample_map)) {
    meta <- sample_map[match(counts$internal_id, sample_map$internal_id), ] %>%
      select(-any_of(c("sample_id", "internal_id", "R1", "R2")))
    out <- bind_cols(out, meta)
  }
  out <- bind_cols(out, tibble(Total_KoRV_reads = total), as_tibble(counts[, cats]))
  list(table = out, rel = as.data.frame(rel, check.names = FALSE), cats = cats,
       n_front = ncol(out) - length(cats))
}

wb <- createWorkbook()
modifyBaseFont(wb, fontSize = 12, fontName = "Arial")
head_style <- createStyle(fontName = "Arial", fontSize = 12, textDecoration = "bold", halign = "center", border = "bottom", borderStyle = "medium")
group_style <- createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold", halign = "center")
center_style <- createStyle(fontName = "Arial", fontSize = 12, halign = "center")

# Column widths that fit the header and the values (capped, so a long sequence column does not take the screen)
fit_widths <- function(df, cap = 60) {
  w <- vapply(seq_along(df), function(j) max(nchar(names(df)[j], type = "width") * 1.1, nchar(format(df[[j]], trim = TRUE), type = "width"), na.rm = TRUE), numeric(1))
  pmin(cap, pmax(8, w + 2))
}
# Arial, centred cells on every sheet; headers bold; widths fitted to the content
finish_sheet <- function(name, nrows, widths) {
  addStyle(wb, name, center_style, rows = seq_len(nrows), cols = seq_along(widths), gridExpand = TRUE, stack = TRUE)
  setColWidths(wb, name, cols = seq_along(widths), widths = widths)
}
plain_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = head_style)
  finish_sheet(name, nrow(df) + 1, fit_widths(df))
}

add_counts_sheet <- function(name, category) {
  s <- summarise_category(category)
  n <- nrow(s$table); k <- length(s$cats); f <- s$n_front
  addWorksheet(wb, name)
  body <- cbind(as.data.frame(s$table, check.names = FALSE), s$rel)
  names(body) <- c(names(s$table), s$cats)
  writeData(wb, name, body, startRow = 4, headerStyle = head_style)
  writeData(wb, name, "Read counts", startCol = f + 1, startRow = 2)
  writeData(wb, name, "Relative abundances", startCol = f + k + 1, startRow = 2)
  mergeCells(wb, name, cols = (f + 1):(f + k), rows = 2)
  mergeCells(wb, name, cols = (f + k + 1):(f + 2 * k), rows = 2)
  addStyle(wb, name, group_style, rows = 2, cols = (f + 1):(f + 2 * k), gridExpand = TRUE)
  foot <- 4 + n + 1
  writeData(wb, name, "Number of samples", startCol = 1, startRow = foot)
  mergeCells(wb, name, cols = 1:f, rows = foot)
  addStyle(wb, name, group_style, rows = foot, cols = 1)
  for (j in seq_len(k)) {
    col <- f + k + j
    L <- int2col(col)
    writeFormula(wb, name, sprintf('COUNTIF(%s5:%s%d,"> 0")', L, L, 4 + n), startCol = col, startRow = foot)
  }
  # white (0) to dark red (100) colour scale on the relative abundances, as in the hand-formatted reports
  conditionalFormatting(wb, name, cols = (f + k + 1):(f + 2 * k), rows = 5:(4 + n), type = "colourScale",
                        style = c("#FFFFFF", "#C00000"), rule = c(0, 100))
  w <- fit_widths(body); w[(f + 1):(f + 2 * k)] <- pmax(w[(f + 1):(f + 2 * k)], 12)
  finish_sheet(name, foot, w)
}

add_counts_sheet(if ("Endogenous" %in% extra_cols) "counts_and_abundances_subtype" else "counts_and_abundances", "Subtype")
if ("Endogenous" %in% extra_cols) add_counts_sheet("counts_and_abundances_sub_endo", "Subtype_Endogenous")

# ---- supporting sheets -------------------------------------------------------------------
fc <- filtered_counts
fc$Sample <- display_name(fc$Sample)
plain_sheet("Reads_filtered_OTU_100", fc)
if (!is.null(stats_file)) {
  st <- as.data.frame(read_tsv(stats_file, col_types = cols(.default = "c"), progress = FALSE), check.names = FALSE)
  st[-1] <- lapply(st[-1], function(x) { y <- suppressWarnings(as.numeric(x)); if (anyNA(y[!is.na(x)])) x else y })   # numbers as numbers
  plain_sheet("sample_read_stats", st)
}
plain_sheet("Representative_KorV_sequences", as.data.frame(assignments, check.names = FALSE))

saveWorkbook(wb, snakemake@output[[1]], overwrite = TRUE)
message("report written for ", length(sample_cols), " samples")
