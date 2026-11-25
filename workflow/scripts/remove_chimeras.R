source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))

remove <- isTRUE(snakemake@params[["drop_chimeras"]])
mode <- snakemake@params[["mode"]]
stopifnot(mode %in% c("legacy", "abundance"))

tab <- read.csv(snakemake@input[["table"]], check.names = FALSE, quote = "", stringsAsFactors = FALSE)
cnt <- read.csv(snakemake@input[["counts"]], check.names = FALSE, quote = "", stringsAsFactors = FALSE)
fa <- readLines(snakemake@input[["fasta"]])
fa_ids <- sub("^>", "", fa[c(TRUE, FALSE)])
fa_seq <- fa[c(FALSE, TRUE)]
chim_lines <- readLines(snakemake@input[["chimeras"]])
chim_ids <- strip_size(sub("^>", "", chim_lines[startsWith(chim_lines, ">")]))
chim_ids <- sub(";size=.*", "", chim_ids)

stopifnot(identical(fa_ids, tab[["#OTU ID"]]))
drop <- if (remove) tab[["#OTU ID"]] %in% chim_ids else rep(FALSE, nrow(tab))

sample_cols <- setdiff(names(tab), c("#OTU ID", "N_samples_present", "Total_number_of_reads"))
chim_reads <- colSums(tab[drop, sample_cols, drop = FALSE])
cnt$Reads_chimeric <- as.integer(chim_reads[match(cnt$Sample, names(chim_reads))])
cnt$Reads_chimeric[is.na(cnt$Reads_chimeric)] <- 0L
cnt$Count_OTU_100_post_filtering <- cnt$Count_OTU_100_post_filtering - cnt$Reads_chimeric

tab <- tab[!drop, , drop = FALSE]
fa_ids <- fa_ids[!drop]
fa_seq <- fa_seq[!drop]

write.csv(tab, snakemake@output[["table"]], quote = FALSE, row.names = FALSE)
write.csv(cnt, snakemake@output[["counts"]], quote = FALSE, row.names = FALSE)
write_fasta_simple(fa_ids, fa_seq, snakemake@output[["fasta"]])
size97 <- if (mode == "abundance") tab$Total_number_of_reads else rep(1L, nrow(tab))
write_fasta_simple(paste0(fa_ids, ";size=", size97), fa_seq, snakemake@output[["fasta97"]])
message(sum(drop), " of ", length(drop), " OTUs removed as chimeras (", sum(chim_reads), " reads); ",
        length(chim_ids), " called chimeric", if (!remove) " (remove_chimeras is false, kept)" else "")
