source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))

# reads -> derep seed (identical sequences)
derep <- read_uc_members(snakemake@input[["derep_uc"]]) %>%
  mutate(sample = sample_of_read(query)) %>%
  count(seed = cluster, sample, name = "count")

# derep seed -> 100% centroid
seed_to_centroid <- read_uc_members(snakemake@input[["cluster_uc"]]) %>%
  rename(seed = query, otu = cluster)

centroids <- Biostrings::readDNAStringSet(snakemake@input[["centroids"]])
otu_ids <- strip_size(names(centroids))
sample_order <- unique(derep$sample)

long <- derep %>%
  inner_join(seed_to_centroid, by = "seed") %>%
  group_by(otu, sample) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(sample_idx = match(sample, sample_order), otu_idx = match(otu, otu_ids)) %>%
  arrange(sample_idx, otu_idx) %>%
  select(otu, sample, count)

stopifnot(!anyNA(long$otu), sum(long$count) == sum(derep$count))
write_tsv(long, snakemake@output[["long"]])
write_tsv(tibble(otu = otu_ids, sequence = toupper(as.character(centroids))),
          snakemake@output[["seqs"]])
message(nrow(long), " OTU x sample cells, ", length(otu_ids), " 100% OTUs, ",
        sum(long$count), " reads")
