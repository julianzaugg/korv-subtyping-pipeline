source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
suppressPackageStartupMessages({library(openxlsx); library(digest)})

otu <- read.csv(snakemake@input[["table"]], check.names = FALSE, quote = "", stringsAsFactors = FALSE)
names(otu)[1] <- "rep_100_ID"
stat_cols <- c("rep_100_ID", "N_samples_present", "Total_number_of_reads")
sample_cols <- setdiff(names(otu), stat_cols)

membership <- read_uc_members(snakemake@input[["uc"]]) %>%
  rename(rep_100_ID = query, rep_97_ID = cluster)
stopifnot(!anyDuplicated(membership$rep_100_ID), all(otu$rep_100_ID %in% membership$rep_100_ID))

# The representative of a 97% cluster is its 100% member with the most reads.
otu <- otu %>%
  left_join(membership, by = "rep_100_ID") %>%
  relocate(rep_97_ID, .after = rep_100_ID) %>%
  arrange(rep_97_ID) %>%
  group_by(rep_97_ID) %>%
  mutate(Cluster_ID = cur_group_id()) %>%
  relocate(Cluster_ID, .after = rep_97_ID) %>%
  ungroup() %>%
  arrange(rep_97_ID, desc(Total_number_of_reads)) %>%
  group_by(rep_97_ID) %>%
  mutate(new_rep_97_ID = rep_100_ID[which.max(Total_number_of_reads)]) %>%
  relocate(new_rep_97_ID, .after = rep_97_ID) %>%
  ungroup() %>%
  as.data.frame(check.names = FALSE)

collapsed <- otu %>%
  group_by(new_rep_97_ID, Cluster_ID) %>%
  summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
  as.data.frame(check.names = FALSE)
collapsed$N_samples_present <- apply(collapsed[, sample_cols, drop = FALSE], 1, function(x) sum(x > 0))
collapsed$Total_number_of_reads <- apply(collapsed[, sample_cols, drop = FALSE], 1, sum)
collapsed <- collapsed %>% relocate(N_samples_present, Total_number_of_reads, .after = Cluster_ID)

write.xlsx(list("OTU_table_100_filt_with_reps" = otu, "OTU_table_100_agg_per_rep_97" = collapsed),
           snakemake@output[["xlsx"]])

seqs <- read_tsv(snakemake@input[["seqs"]], col_types = "cc", progress = FALSE)
final <- collapsed %>%
  select(new_rep_97_ID, Cluster_ID, N_samples_present, Total_number_of_reads) %>%
  left_join(seqs, by = c("new_rep_97_ID" = "otu")) %>%
  mutate(seq_md5 = vapply(sequence, digest, "", algo = "md5", serialize = FALSE))
write_fasta_simple(final$new_rep_97_ID, final$sequence, snakemake@output[["fasta"]])
write_tsv(final, snakemake@output[["tsv"]])
message(nrow(final), " final 97% representatives from ", nrow(otu), " 100% OTUs")
