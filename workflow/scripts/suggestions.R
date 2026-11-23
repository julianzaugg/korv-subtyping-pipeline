source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
source(file.path(dirname(snakemake@scriptdir), "scripts", "subtype_lib.R"))
suppressPackageStartupMessages({library(openxlsx); library(yaml); library(digest); library(stringr)})

rules <- read_yaml(snakemake@input[["rules"]])
refs <- readAAStringSet(snakemake@input[["refs"]])
ref_names <- names(refs); ref_safe <- safe_name(ref_names); ref_sub <- reference_subtype(ref_names)
final <- read_tsv(snakemake@input[["final"]], col_types = cols(.default = "c"), progress = FALSE) %>%
  mutate(Total_number_of_reads = as.numeric(Total_number_of_reads), N_samples_present = as.numeric(N_samples_present))
tr <- read_tsv(snakemake@input[["translation"]], col_types = cols(.default = "c"), progress = FALSE) %>%
  mutate(n_internal_stops = as.integer(n_internal_stops))
lib <- read_tsv(snakemake@input[["library"]], col_types = cols(.default = "c"), progress = FALSE)
hits <- tryCatch(read_tsv(snakemake@input[["hits"]], col_names = c("query", "md5", "id", "qcov", "tcov"),
                          col_types = "ccddd", progress = FALSE), error = function(e) tibble(query = character(), md5 = character(), id = numeric(), qcov = numeric(), tcov = numeric()))
chim_ids <- sub(";size=.*", "", names(readDNAStringSet(snakemake@input[["chimeras"]])))

# ---- alignment ------------------------------------------------------------------------------
aln <- toupper(as.matrix(readAAStringSet(snakemake@input[["aln"]])))
ref_rows <- match(ref_safe, rownames(aln)); stopifnot(!anyNA(ref_rows))
R <- aln[ref_rows, , drop = FALSE]
ref_present <- colMeans(R != "-")                      # fraction of references with a residue per column

nq <- anchored_region_columns(R, ref_safe, rules$nqqhnn_zone)
vr <- anchored_region_columns(R, ref_safe, rules$variable_region)
if (length(nq) == 0) warning("NQQHNN anchors not found; the zone is treated as conserved")
if (length(vr) == 0) warning("variable region anchors not found; it is treated as conserved")
in_nq <- seq_len(ncol(aln)) %in% nq
in_vr <- seq_len(ncol(aln)) %in% vr

# ---- per-sequence evidence --------------------------------------------------------------------
ev_rows <- lapply(final$new_rep_97_ID, function(rep_id) {
  x <- aln[rep_id, ]; ng <- x != "-"
  first <- min(which(ng)); last <- max(which(ng))
  both <- sweep(R != "-", 2, ng, "&"); match <- both & sweep(R, 2, x, "==")
  id_all <- rowSums(match) / pmax(rowSums(both), 1)
  id_nq <- if (length(nq)) rowSums(match[, nq, drop = FALSE]) / pmax(rowSums(both[, nq, drop = FALSE]), 1) else rep(NA_real_, nrow(R))
  id_vr <- if (length(vr)) rowSums(match[, vr, drop = FALSE]) / pmax(rowSums(both[, vr, drop = FALSE]), 1) else rep(NA_real_, nrow(R))
  cons <- !(in_nq | in_vr)    # conserved columns: outside the NQQHNN zone and the variable region
  id_cons <- rowSums(match[, cons, drop = FALSE]) / pmax(rowSums(both[, cons, drop = FALSE]), 1)
  cover <- rowSums(both) / rowSums(R != "-")
  best <- which.max(id_all)
  other <- which(ref_sub != ref_sub[best]); best_other <- other[which.max(id_all[other])]
  interior <- seq_along(x) >= first & seq_along(x) <= last
  # missing / inserted sequence is measured against this sequence's own nearest reference, so a
  # sequence that matches its subtype exactly is never penalised for features of other subtypes
  nr <- R[best, ] != "-"
  del_cols <- !ng & interior & nr
  tibble(new_rep_97_ID = rep_id,
         nearest_reference = ref_names[best], nearest_reference_subtype = ref_sub[best], nearest_reference_identity = id_all[best],
         nearest_reference_coverage = cover[best],
         other_subtype_reference = ref_names[best_other], other_subtype_identity = id_all[best_other],
         margin_to_other_subtype = id_all[best] - id_all[best_other],
         conserved_identity = max(id_cons),
         nqqhnn_nearest_reference = if (length(nq)) ref_names[which.max(id_nq)] else NA_character_,
         nqqhnn_identity = if (length(nq)) max(id_nq) else NA_real_,
         variable_region_nearest_reference = if (length(vr)) ref_names[which.max(id_vr)] else NA_character_,
         variable_region_identity = if (length(vr)) max(id_vr) else NA_real_,
         terminal_missing_5_bp = 3 * sum(seq_along(x) < first & nr),
         terminal_missing_3_bp = 3 * sum(seq_along(x) > last & nr),
         deletion_conserved_bp = 3 * sum(del_cols & !in_nq & !in_vr),
         missing_nqqhnn_bp = 3 * sum(del_cols & in_nq),
         variable_region_missing_bp = 3 * sum(del_cols & in_vr),    # informational
         variable_region_residues = sum(ng & in_vr),                # flagged only if very short, see rules
         inserted_residues = sum(ng & !nr),
         across = list(setNames(id_all, ref_names)))
})
ev <- bind_rows(ev_rows)
ident_mat <- do.call(rbind, ev$across); rownames(ident_mat) <- ev$new_rep_97_ID
ev$across <- NULL

# ---- library evidence -------------------------------------------------------------------------
sg <- rules$suggestion
lib_labels <- lib %>% group_by(seq_md5) %>%
  summarise(labels = paste(unique(Subtype), collapse = " | "), endo = paste(na.omit(unique(Endogenous[Endogenous != ""])), collapse = "|"),
            sources = paste(unique(source), collapse = ","), .groups = "drop")
md5_oriented <- setNames(vapply(tr$nt_oriented, digest, "", algo = "md5", serialize = FALSE), tr$id)
hits <- hits %>% mutate(id = id / 100) %>% inner_join(lib_labels, by = c("md5" = "seq_md5")) %>% arrange(query, desc(id))
lib_ev <- lapply(final$new_rep_97_ID, function(rep_id) {
  h <- hits %>% filter(query == rep_id)
  exact <- md5_oriented[[rep_id]] %in% lib_labels$seq_md5
  if (exact) h <- bind_rows(tibble(query = rep_id, md5 = md5_oriented[[rep_id]], id = 1, qcov = 1, tcov = 1) %>% inner_join(lib_labels, by = c("md5" = "seq_md5")), h %>% filter(md5 != md5_oriented[[rep_id]]))
  if (nrow(h) == 0) return(tibble(new_rep_97_ID = rep_id, library_exact_match = FALSE, library_best_identity = NA_real_, library_labels_best = NA_character_,
                                  library_labels_window = NA_character_, library_n_window_labels = 0L, library_majority_label = NA_character_,
                                  library_top_hits = NA_character_, library_endogenous = NA_character_))
  w <- h %>% filter(id >= h$id[1] - sg$label_agreement_window)
  labs <- unlist(strsplit(w$labels, " \\| ", perl = TRUE))
  tab <- sort(table(labs), decreasing = TRUE)
  top <- head(h, sg$top_hits)
  tibble(new_rep_97_ID = rep_id, library_exact_match = exact, library_best_identity = h$id[1], library_labels_best = h$labels[1],
         library_labels_window = paste(unique(labs), collapse = " | "), library_n_window_labels = length(unique(labs)),
         library_majority_label = names(tab)[1],
         library_top_hits = paste(sprintf("%s (%.1f%%; %s)", top$labels, 100 * top$id, top$sources), collapse = "; "),
         library_endogenous = h$endo[1])
})
ev <- ev %>% left_join(bind_rows(lib_ev), by = "new_rep_97_ID")

# ---- flags and suggestion ---------------------------------------------------------------------
ev <- ev %>% left_join(tr %>% select(new_rep_97_ID = id, strand, frame, n_internal_stops, protein), by = "new_rep_97_ID") %>%
  mutate(chimera = new_rep_97_ID %in% chim_ids)
ev$flags <- apply(ev, 1, function(r) {
  f <- character()
  if (as.integer(r[["n_internal_stops"]]) > rules$internal_stops$max) f <- c(f, "internal_stop")
  if (as.numeric(r[["deletion_conserved_bp"]]) > rules$conserved_regions$max_deletion_bp) f <- c(f, "deletion_conserved_region")
  if (max(as.numeric(r[["terminal_missing_5_bp"]]), as.numeric(r[["terminal_missing_3_bp"]])) > rules$conserved_regions$max_terminal_missing_bp) f <- c(f, "end_missing")
  if (as.numeric(r[["missing_nqqhnn_bp"]]) > rules$nqqhnn_zone$max_missing_bp) f <- c(f, "nqqhnn_deletion")
  if (length(vr) && !is.null(rules$variable_region$min_residues) && as.numeric(r[["variable_region_residues"]]) < rules$variable_region$min_residues) f <- c(f, "variable_region_short")
  if (r[["chimera"]] == "TRUE") f <- c(f, "possible_chimera")
  paste(f, collapse = ", ")
})
nf_flags <- c("internal_stop", "deletion_conserved_region", "end_missing", "nqqhnn_deletion", "variable_region_short")
ev$nonfunctional_features <- vapply(strsplit(ev$flags, ", "), function(f) paste(intersect(f, nf_flags), collapse = ", "), "")

decide <- function(r) {
  nf <- r[["nonfunctional_features"]] != ""
  lib_id <- suppressWarnings(as.numeric(r[["library_best_identity"]])); ref_id <- as.numeric(r[["nearest_reference_identity"]])
  nlab <- as.integer(r[["library_n_window_labels"]]); lab <- r[["library_majority_label"]]
  lib_ok <- !is.na(lib_id) && lib_id >= sg$check_library_identity
  cons_id <- as.numeric(r[["conserved_identity"]])
  novel <- (is.na(lib_id) || lib_id < rules$novelty$min_library_identity) && cons_id < rules$novelty$min_conserved_identity
  if (lib_ok) {
    if (nlab == 1) {
      tier <- if (lib_id >= sg$high_library_identity) "high" else "check"
      basis <- sprintf("curated library, %.1f%% nucleotide identity", 100 * lib_id)
    } else {
      tier <- "check"; basis <- sprintf("library hits disagree (%s), majority label", r[["library_labels_window"]])
    }
    if (nf && !grepl("non-?functional", lab, ignore.case = TRUE) && tier == "high") { tier <- "check"; basis <- paste0(basis, "; but non-functional features present") }
    return(c(lab, tier, basis))
  }
  # a sequence with non-functional features is not a functional new subtype, whatever its distance to the references
  if (nf) return(c("Non-functional", "check", paste0("non-functional criteria: ", r[["nonfunctional_features"]],
                   if (novel) sprintf(" (also far from every reference, %.1f%% protein identity in the conserved regions)", 100 * cons_id) else "")))
  if (novel) return(c(sprintf("possible new subtype (nearest %s)", r[["nearest_reference_subtype"]]), "novel",
                      sprintf("%.1f%% protein identity to nearest reference in the conserved regions, no close library sequence", 100 * cons_id)))
  c(r[["nearest_reference_subtype"]], "check", sprintf("nearest A-N reference (%.1f%% protein identity), no close library sequence", 100 * ref_id))
}
dec <- t(apply(ev, 1, decide))
ev$suggested_subtype <- dec[, 1]; ev$confidence <- dec[, 2]; ev$basis <- dec[, 3]

# ---- output ---------------------------------------------------------------------------------------
collapsed <- read.xlsx(snakemake@input[["collapsed"]], sheet = "OTU_table_100_agg_per_rep_97", check.names = FALSE)
samp <- setdiff(names(collapsed), c("new_rep_97_ID", "Cluster_ID", "N_samples_present", "Total_number_of_reads"))
top_sample <- samp[apply(as.matrix(collapsed[, samp, drop = FALSE]), 1, which.max)]
top <- tibble(new_rep_97_ID = collapsed$new_rep_97_ID, top_sample = top_sample)

out <- final %>% select(new_rep_97_ID, seq_md5, Total_number_of_reads, N_samples_present, nt_sequence = sequence) %>%
  left_join(top, by = "new_rep_97_ID") %>% left_join(ev, by = "new_rep_97_ID") %>%
  mutate(Subtype = "", Endogenous = "", Unsure = "", Note = "") %>%
  arrange(factor(confidence, c("novel", "check", "high")), suggested_subtype, desc(Total_number_of_reads))
sheet <- out %>% transmute(
  sample_id = top_sample, new_rep_97_ID, seq_md5, Total_number_of_reads, N_samples_present,
  Subtype, Endogenous, Unsure, Note,
  suggested_subtype, confidence, basis, flags,
  nearest_reference, nearest_reference_identity = round(nearest_reference_identity, 4),
  other_subtype_reference, other_subtype_identity = round(other_subtype_identity, 4), margin_to_other_subtype = round(margin_to_other_subtype, 4),
  conserved_identity = round(conserved_identity, 4),
  variable_region_nearest_reference, variable_region_identity = round(variable_region_identity, 4),
  nqqhnn_nearest_reference, nqqhnn_identity = round(nqqhnn_identity, 4),
  library_exact_match, library_best_identity = round(library_best_identity, 4), library_labels_best, library_top_hits, library_endogenous,
  n_internal_stops, terminal_missing_5_bp, terminal_missing_3_bp, deletion_conserved_bp, missing_nqqhnn_bp, variable_region_missing_bp, variable_region_residues, inserted_residues,
  strand, frame, Sequence = sub("\\*$", "", protein))

legend <- tibble(column = c("Subtype, Endogenous, Unsure, Note", "suggested_subtype / confidence / basis", "flags",
                            "nearest_reference*", "other_subtype_reference*", "variable_region_*", "nqqhnn_*", "library_*", "terminal_missing / deletion_conserved / missing_nqqhnn / inserted", "sample_id"),
  meaning = c("Yellow columns, fill these in. Only Subtype is required. Endogenous is optional free text (blank column is ignored, blanks in a partly filled column become unassigned) and gets its own summary sheet. Unsure and Note are kept but not used. See docs/outputs.md.",
              "A suggestion only. high: close match to curated library sequences that agree. check: needs a look. novel: rare, below 80% identity to every reference in the conserved regions and no close library sequence.",
              "internal_stop, deletion_conserved_region, end_missing, nqqhnn_deletion, variable_region_short (non-functional criteria from resources/nonfunctional_rules.yaml), possible_chimera. Only a very short variable region is flagged, other variable region indels are not.",
              "Nearest A-N reference by protein identity over the columns where both sequences have a residue, and that identity. conserved_identity is the best identity over the conserved columns only (outside the NQQHNN zone and the variable region).",
              "Best reference of a different subtype, and the identity margin to it (small margin = ambiguous).",
              "variable_region_*: nearest reference and identity within the variable region (between ASQQA and LGCSY in the A reference) only; it needs human review. Only a very short region (variable_region_residues below the minimum) is flagged.",
              "Nearest reference and identity within the NQQHNN zone near the start (between GGAEM and PHQPM in the A reference). This zone is variable and some missing residues are allowed.",
              "Nearest curated sequences from resources/subtype_library.tsv (nucleotide identity, label, source analysis).",
              "Missing or inserted sequence in base pairs / residues relative to the references.",
              "Sample with most reads for this sequence."))
wb <- createWorkbook()
addWorksheet(wb, "suggestions"); writeData(wb, "suggestions", sheet, headerStyle = createStyle(textDecoration = "bold", border = "bottom"))
freezePane(wb, "suggestions", firstActiveRow = 2, firstActiveCol = 4)
addStyle(wb, "suggestions", createStyle(fgFill = "#FFF2CC"), rows = 2:(nrow(sheet) + 1), cols = 6:9, gridExpand = TRUE, stack = TRUE)
conditionalFormatting(wb, "suggestions", cols = 11, rows = 2:(nrow(sheet) + 1), rule = "novel", type = "contains", style = createStyle(bgFill = "#F4B6B6"))
conditionalFormatting(wb, "suggestions", cols = 11, rows = 2:(nrow(sheet) + 1), rule = "check", type = "contains", style = createStyle(bgFill = "#FCE4B6"))
conditionalFormatting(wb, "suggestions", cols = 11, rows = 2:(nrow(sheet) + 1), rule = "high", type = "contains", style = createStyle(bgFill = "#CFE8CF"))
setColWidths(wb, "suggestions", cols = 1:ncol(sheet), widths = "auto")
addWorksheet(wb, "legend"); writeData(wb, "legend", legend)
setColWidths(wb, "legend", cols = 1:2, widths = c(45, 120))
saveWorkbook(wb, snakemake@output[["xlsx"]], overwrite = TRUE)

ident_df <- as_tibble(ident_mat, rownames = "new_rep_97_ID")
names(ident_df)[-1] <- paste0("identity__", names(ident_df)[-1])
write_tsv(out %>% select(-nt_sequence, -Subtype, -Endogenous, -Unsure, -Note) %>% left_join(ident_df, by = "new_rep_97_ID"), snakemake@output[["evidence"]])
message(nrow(out), " sequences; confidence: ", paste(names(table(out$confidence)), table(out$confidence), collapse = ", "))
