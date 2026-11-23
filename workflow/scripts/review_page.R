# Interactive review page: one self-contained HTML file with the tree and the whole alignment side
# by side, linked selection, and subtype assignment that exports the completed suggestions workbook.
# The data are embedded as JSON; the page (workflow/scripts/review_page/) needs no server or network.
source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
source(file.path(dirname(snakemake@scriptdir), "scripts", "subtype_lib.R"))
suppressPackageStartupMessages({library(openxlsx); library(ape); library(yaml); library(jsonlite); library(digest)})

rules <- read_yaml(snakemake@input[["rules"]])
sheet <- read.xlsx(snakemake@input[["suggestions"]], sheet = "suggestions", check.names = FALSE)
legend <- read.xlsx(snakemake@input[["suggestions"]], sheet = "legend", check.names = FALSE)
ev <- read_tsv(snakemake@input[["evidence"]], col_types = cols(.default = "c"), progress = FALSE)
collapsed <- read.xlsx(snakemake@input[["collapsed"]], sheet = "OTU_table_100_agg_per_rep_97", check.names = FALSE)
analysis_name <- snakemake@params[["analysis_name"]]

ref_full <- sub("^identity__", "", grep("^identity__", names(ev), value = TRUE))
ref_safe <- safe_name(ref_full); ref_sub <- reference_subtype(ref_full)

# ---- alignment --------------------------------------------------------------------------------------
aln <- toupper(as.matrix(Biostrings::readAAStringSet(snakemake@input[["aln"]])))
stopifnot(all(sheet$new_rep_97_ID %in% rownames(aln)), all(ref_safe %in% rownames(aln)))
row_ids <- rownames(aln)
Rm <- aln[match(ref_safe, rownames(aln)), , drop = FALSE]
vr_cols <- anchored_region_columns(Rm, ref_safe, rules$variable_region)
nq_cols <- anchored_region_columns(Rm, ref_safe, rules$nqqhnn_zone)
fills <- clustalx_fills(aln)
codes <- c("#80a0f0" = "b", "#f01505" = "r", "#c048c0" = "m", "#15c015" = "g", "#f09048" = "o", "#15a4a4" = "c", "#f08080" = "p", "#c0c000" = "y", "#ffffff" = ".")
fill_codes <- apply(fills, 1, function(x) paste(codes[x], collapse = ""))
conservation <- apply(aln, 2, function(x) { r <- x[x != "-"]; if (length(r) == 0) 0 else max(table(r)) / length(r) })

# ---- tree (edge table, so the page needs no Newick parser) ---------------------------------------------
tree <- ape::read.tree(snakemake@input[["treefile"]])
stopifnot(setequal(tree$tip.label, row_ids))
ntip <- length(tree$tip.label); nn <- ntip + tree$Nnode
parent <- rep(-1L, nn); bl <- rep(0, nn)
parent[tree$edge[, 2]] <- tree$edge[, 1] - 1L
bl[tree$edge[, 2]] <- if (is.null(tree$edge.length)) 0 else tree$edge.length
nodes <- lapply(seq_len(nn), function(i) list(parent = parent[i], bl = bl[i], tip = if (i <= ntip) tree$tip.label[i] else NA_character_))

# ---- sequences ---------------------------------------------------------------------------------------
ev_by_id <- split(seq_len(nrow(ev)), ev$new_rep_97_ID)
sample_cols <- setdiff(names(collapsed), c("new_rep_97_ID", "Cluster_ID", "N_samples_present", "Total_number_of_reads"))
num <- function(x) suppressWarnings(as.numeric(x))
sequences <- lapply(seq_len(nrow(sheet)), function(i) {
  id <- sheet$new_rep_97_ID[i]; e <- ev[ev_by_id[[id]][1], ]
  cr <- collapsed[match(id, collapsed$new_rep_97_ID), sample_cols, drop = FALSE]
  per_sample <- as.list(unlist(cr)[unlist(cr) > 0])
  flags <- if (is.na(sheet$flags[i]) || sheet$flags[i] == "") character() else strsplit(sheet$flags[i], ", ", fixed = TRUE)[[1]]
  list(id = id, group = group_of(sheet$suggested_subtype[i]), nearest_ref_row = safe_name(sheet$nearest_reference[i]),
       flags = I(flags), per_sample = per_sample,
       identity = as.list(setNames(num(unlist(e[paste0("identity__", ref_full)])), ref_full)),
       export = as.list(sheet[i, , drop = FALSE]))
})

groups <- sort(unique(c(vapply(sequences, `[[`, "", "group"), ref_sub)))
pal <- subtype_palette(groups)
vocab <- unique(c(sort(unique(ref_sub)), "Non-functional"))
fingerprint <- digest(paste(sort(paste(sheet$new_rep_97_ID, sheet$seq_md5)), collapse = "|"), algo = "md5", serialize = FALSE)

data <- list(
  schema_version = 1L, analysis_name = analysis_name, generated = format(Sys.Date()), fingerprint = fingerprint,
  rules = list(variable_region_cols = I(as.integer(vr_cols)), nqqhnn_cols = I(as.integer(nq_cols)),
               variable_region_min_residues = rules$variable_region$min_residues),
  palette = as.list(pal),
  subtype_order = if (is.null(rules$subtype_order)) NULL else list(
    order = I(as.character(rules$subtype_order$order)), intermediates = rules$subtype_order$intermediates),
  references = lapply(seq_along(ref_full), function(i) list(id = ref_safe[i], name = ref_full[i], subtype = ref_sub[i])),
  tree = list(nodes = nodes),
  alignment = list(ncol = ncol(aln), order = I(row_ids),
                   seq = as.list(setNames(apply(aln, 1, paste, collapse = ""), row_ids)),
                   fill = as.list(setNames(fill_codes[row_ids], row_ids)),
                   conservation = I(round(conservation, 3))),
  sequences = sequences,
  export_columns = I(names(sheet)),
  legend = lapply(seq_len(nrow(legend)), function(i) list(column = legend[[1]][i], meaning = legend[[2]][i])),
  vocabulary_seed = I(vocab))

json <- toJSON(data, auto_unbox = TRUE, na = "null", null = "null", digits = NA)
json <- gsub("</", "<\\u002f", json, fixed = TRUE)          # keeps the JSON from closing the script tag

dir_t <- file.path(dirname(snakemake@scriptdir), "scripts", "review_page")
rd <- function(f) paste(readLines(file.path(dir_t, f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
fill <- function(tpl, key, val) { p <- regexpr(key, tpl, fixed = TRUE); stopifnot(p > 0); paste0(substr(tpl, 1, p - 1), val, substr(tpl, p + nchar(key), nchar(tpl))) }
html <- rd("template.html")
html <- fill(html, "__TITLE__", paste("KoRV subtype review:", analysis_name))
html <- fill(html, "__STYLE__", rd("app.css"))
html <- fill(html, "__CORE__", rd("core.js"))
html <- fill(html, "__APP__", rd("app.js"))
html <- fill(html, "__DATA__", json)           # last: the data cannot contain a placeholder that matters
writeLines(html, snakemake@output[["html"]], useBytes = TRUE)
message("review page written: ", nrow(sheet), " sequences, ", length(row_ids), " rows, ", round(file.size(snakemake@output[["html"]]) / 1024), " KB")
