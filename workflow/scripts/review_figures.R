source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
source(file.path(dirname(snakemake@scriptdir), "scripts", "subtype_lib.R"))
suppressPackageStartupMessages({library(ggplot2); library(ggtree); library(ape); library(yaml); library(stringr)})

rules <- read_yaml(snakemake@input[["rules"]])
ev <- read_tsv(snakemake@input[["evidence"]], col_types = cols(.default = "c"), progress = FALSE) %>%
  mutate(Total_number_of_reads = as.numeric(Total_number_of_reads), N_samples_present = as.numeric(N_samples_present))
dir.create(dirname(snakemake@output[["tree_cladogram"]]), showWarnings = FALSE, recursive = TRUE)

# PDF with base pdf() and PNG with ragg when it is installed: cairo devices space letters unevenly.
save_fig <- function(plot, pdf_path, width, height, dpi = 110) {
  ggsave(pdf_path, plot, width = width, height = height, limitsize = FALSE,
         device = function(filename, ...) grDevices::pdf(filename, ..., useDingbats = FALSE))
  dev <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(sub("\\.pdf$", ".png", pdf_path), plot, width = width, height = height, dpi = dpi,
         limitsize = FALSE, device = dev)
}

ev$group <- group_of(ev$suggested_subtype)
ev$label <- sprintf("%s  [%s, %s reads, %s samples]", ev$new_rep_97_ID, ev$group, ev$Total_number_of_reads, ev$N_samples_present)
ev$flagged <- ev$flags != "" & !is.na(ev$flags)

ref_full <- sub("^identity__", "", grep("^identity__", names(ev), value = TRUE))
ref_safe <- safe_name(ref_full); ref_sub <- reference_subtype(ref_full)

groups <- sort(unique(c(ev$group, ref_sub)))
pal <- subtype_palette(groups)

# ---- tree -------------------------------------------------------------------------------------------
tree <- ape::read.tree(snakemake@input[["treefile"]])
meta <- tibble(tip = tree$tip.label) %>%
  left_join(ev %>% select(tip = new_rep_97_ID, group, tip_label = label, flagged), by = "tip") %>%
  mutate(is_ref = tip %in% ref_safe)
refmeta <- tibble(tip = ref_safe, rgroup = ref_sub, rlabel = ref_full)
meta <- meta %>% left_join(refmeta, by = "tip") %>%
  mutate(group = ifelse(is_ref, rgroup, group), tip_label = ifelse(is_ref, paste0("REF ", rlabel), tip_label), flagged = ifelse(is_ref, FALSE, flagged))
# Two views of the same tree: a cladogram (topology only, equal-length branches) and a
# distance tree (branch lengths are expected amino acid substitutions per site).
# A triangle marks a sequence with at least one flag (see the flags column of the suggestions workbook).
n <- length(tree$tip.label)
save_tree <- function(branch_length, path) {
  g0 <- ggtree(tree, branch.length = branch_length)
  mx <- max(g0$data$x)
  p <- g0 %<+% meta +
    geom_tiplab(aes(label = tip_label, color = group, fontface = ifelse(is_ref, "bold", "plain")), size = 2.4, show.legend = FALSE) +
    geom_tippoint(aes(color = group, shape = flagged), size = 2) +
    scale_color_manual(values = pal, name = "subtype (reference or suggested)", na.translate = FALSE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "flags raised", labels = c("none", "one or more"), na.translate = FALSE) +
    xlim(0, mx * 2.4) + theme(legend.position = "left")
  if (branch_length != "none") {
    p <- p + geom_treescale(fontsize = 2.8) +
      labs(caption = "Scale bar: expected amino acid substitutions per site")
  }
  save_fig(p, path, width = 14, height = max(6, 0.14 * n + 2))
}
save_tree("none", snakemake@output[["tree_cladogram"]])
save_tree("branch.length", snakemake@output[["tree_distance"]])

# ---- identity heatmap -----------------------------------------------------------------------------------
idm <- as.matrix(ev[, paste0("identity__", ref_full)]); mode(idm) <- "numeric"; rownames(idm) <- ev$new_rep_97_ID
ord <- if (nrow(idm) > 2) hclust(dist(idm))$order else seq_len(nrow(idm))
long <- tibble(rep = rep(ev$new_rep_97_ID, times = length(ref_full)), ref = rep(ref_full, each = nrow(ev)), identity = as.vector(idm)) %>%
  left_join(ev %>% select(rep = new_rep_97_ID, label), by = "rep") %>%
  mutate(label = factor(label, levels = ev$label[ord]), ref = factor(ref, levels = ref_full[order(ref_sub)]))
h <- ggplot(long, aes(ref, label, fill = identity)) + geom_tile(color = "white", linewidth = 0.1) +
  geom_text(aes(label = sprintf("%.0f", 100 * identity), color = identity > 0.88), size = 1.8, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0.7, 1), oob = scales::squish, name = "Protein identity") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Protein identity of each sequence to the A-N references") +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
# width follows the number of reference columns (about 0.17 in each) so cells stay compact
hm_width <- 4.8 + 0.17 * length(ref_full)
save_fig(h, snakemake@output[["heatmap"]], width = hm_width, height = max(4, 0.13 * nrow(ev) + 2))

# ---- variable-region alignment --------------------------------------------------------------------------
aln <- toupper(as.matrix(Biostrings::readAAStringSet(snakemake@input[["aln"]])))
Rm <- aln[match(ref_safe, rownames(aln)), , drop = FALSE]
region <- anchored_region_columns(Rm, ref_safe, rules$variable_region)
if (length(region) > 0) {
  cols <- max(1, min(region) - 4):min(ncol(aln), max(region) + 4)
  # rows follow the order of the alignment file: mafft --reorder writes sequences in guide-tree
  # order, so similar sequences (references and representatives alike) are neighbours
  rows <- rownames(aln)
  nearest <- setNames(safe_name(ev$nearest_reference), ev$new_rep_97_ID)
  df <- do.call(rbind, lapply(rows, function(rn) {
    x <- aln[rn, cols]; is_ref <- rn %in% ref_safe
    ref_x <- if (is_ref) x else aln[nearest[[rn]], cols]
    state <- if (is_ref) rep("reference", length(x)) else ifelse(x == "-", "gap", ifelse(x == ref_x, "same as nearest reference", "differs"))
    tibble(row = rn, col = seq_along(cols), residue = x, state = state, in_region = cols %in% region)
  }))
  lab <- c(setNames(paste0("REF ", ref_full), ref_safe), setNames(ev$label, ev$new_rep_97_ID))
  stopifnot(setequal(rows, names(lab)))
  df$row <- factor(lab[df$row], levels = rev(lab[rows]))
  g <- ggplot(df, aes(col, row, fill = state)) + geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = residue), size = 2.2) +
    scale_fill_manual(values = c(reference = "#cfe3f5", `same as nearest reference` = "#f2f2f2", differs = "#f6b26b", gap = "#ffffff"), name = NULL) +
    annotate("rect", xmin = min(which(cols %in% region)) - 0.5, xmax = max(which(cols %in% region)) + 0.5, ymin = 0.5, ymax = length(rows) + 0.5,
             fill = NA, color = "black", linewidth = 0.5) +
    scale_x_continuous(expand = c(0, 0)) + guides(fill = guide_legend(override.aes = list(color = "grey60"))) +
    labs(x = "boxed columns: variable region", y = NULL, title = "Variable region cut from the full alignment, rows in alignment order") +
    theme_minimal(base_size = 8) + theme(panel.grid = element_blank(), axis.text.x = element_blank())
  save_fig(g, snakemake@output[["variable"]], width = 6 + 0.22 * length(cols), height = max(4, 0.14 * length(rows) + 2))
} else {
  pdf(snakemake@output[["variable"]]); plot.new(); text(0.5, 0.5, "variable region anchors not found in alignment"); dev.off()
}

# ---- whole alignment, ClustalX colours -------------------------------------------------------------------
# One image with every sequence and every column, in alignment order (no paging, however many
# sequences there are). The top row shades each column by the share of its most common residue.
rows <- rownames(aln)
lab_all <- c(setNames(paste0("REF ", ref_full), ref_safe), setNames(ev$label, ev$new_rep_97_ID))
stopifnot(setequal(rows, names(lab_all)))
fills <- clustalx_fills(aln)
conservation <- apply(aln, 2, function(x) { r <- x[x != "-"]; if (length(r) == 0) 0 else max(table(r)) / length(r) })
cons_fill <- colorRampPalette(c("#ffffff", "#303030"))(101)[round(conservation * 100) + 1]
row_levels <- c(rev(lab_all[rows]), "Conservation")           # first level is drawn at the bottom
adf <- rbind(
  tibble(row = rep(rows, times = ncol(aln)), col = rep(seq_len(ncol(aln)), each = length(rows)),
         residue = as.vector(aln[rows, ]), fill = as.vector(fills[rows, ])) %>% mutate(row = lab_all[row]),
  tibble(row = "Conservation", col = seq_len(ncol(aln)), residue = "", fill = cons_fill))
adf$row <- factor(adf$row, levels = row_levels)
is_ref_row <- c(rev(rows %in% ref_safe), FALSE)
ag <- ggplot(adf, aes(col, row, fill = fill)) + geom_tile(color = "white", linewidth = 0.1) +
  geom_text(aes(label = residue), size = 2.3) +
  scale_fill_identity() +
  scale_x_continuous(breaks = seq(10, ncol(aln), by = 10), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = parse(text = sprintf('%s("%s")', ifelse(is_ref_row, "bold", "plain"), row_levels))) +
  labs(x = "alignment column", y = NULL, title = "Alignment of representatives and references",
       caption = "ClustalX colours: blue hydrophobic, red positive, magenta negative, green polar, pink cysteine, orange glycine, yellow proline, cyan aromatic. Top row: share of the most common residue.") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank())
if (length(region) > 0) {
  ag <- ag + annotate("rect", xmin = min(region) - 0.5, xmax = max(region) + 0.5, ymin = 0.5, ymax = length(row_levels) + 0.5,
                      fill = NA, color = "black", linewidth = 0.5)
}
save_fig(ag, snakemake@output[["alignment"]], width = 4.5 + 0.105 * ncol(aln), height = 0.14 * length(row_levels) + 1.2)
message("figures written")
