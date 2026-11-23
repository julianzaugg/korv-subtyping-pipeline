source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))

min_reads <- as.integer(snakemake@params[["min_reads"]])
mode <- snakemake@params[["mode"]]
stopifnot(mode %in% c("legacy", "abundance"))

long <- read_tsv(snakemake@input[["long"]], col_types = "ccd", progress = FALSE)
seqs <- read_tsv(snakemake@input[["seqs"]], col_types = "cc", progress = FALSE)
sample_order <- unique(long$sample)

totals <- long %>%
  group_by(otu) %>%
  summarise(N_samples_present = sum(count > 0), Total_number_of_reads = sum(count), .groups = "drop")
keep <- totals %>% filter(Total_number_of_reads >= min_reads)

wide <- long %>%
  semi_join(keep, by = "otu") %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0) %>%
  left_join(totals, by = "otu") %>%
  mutate(.i = match(otu, seqs$otu)) %>%
  arrange(.i) %>%
  select(`#OTU ID` = otu, N_samples_present, Total_number_of_reads, all_of(sample_order)) %>%
  as.data.frame(check.names = FALSE)

write.csv(wide, snakemake@output[["table"]], quote = FALSE, row.names = FALSE)

per_sample <- tibble(Sample = sample_order) %>%
  left_join(long %>% group_by(Sample = sample) %>% summarise(Count_OTU_100_pre_filtering = sum(count)),
            by = "Sample") %>%
  left_join(long %>% semi_join(keep, by = "otu") %>% group_by(Sample = sample) %>%
              summarise(Count_OTU_100_post_filtering = sum(count)), by = "Sample") %>%
  mutate(Count_OTU_100_post_filtering = coalesce(Count_OTU_100_post_filtering, 0),
         Reads_filtered = Count_OTU_100_pre_filtering - Count_OTU_100_post_filtering)
write.csv(per_sample, snakemake@output[["counts"]], quote = FALSE, row.names = FALSE)

kept_seqs <- seqs[match(wide[["#OTU ID"]], seqs$otu), ]
write_fasta_simple(kept_seqs$otu, kept_seqs$sequence, snakemake@output[["fasta"]])

size97 <- if (mode == "abundance") wide$Total_number_of_reads else rep(1L, nrow(wide))
write_fasta_simple(paste0(kept_seqs$otu, ";size=", size97), kept_seqs$sequence,
                   snakemake@output[["fasta97"]])
message(nrow(wide), " of ", nrow(totals), " OTUs kept (>= ", min_reads, " reads)")
