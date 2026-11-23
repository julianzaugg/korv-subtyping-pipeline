source(file.path(dirname(snakemake@scriptdir), "scripts", "common.R"))
source(file.path(dirname(snakemake@scriptdir), "scripts", "subtype_lib.R"))

nt <- readDNAStringSet(snakemake@input[["fasta"]])
names(nt) <- sub(" .*", "", names(nt))
refs <- readAAStringSet(snakemake@input[["refs"]])

tr <- translate_best_frame(nt, refs, cores = max(1, snakemake@threads)) %>%
  mutate(n_internal_stops = n_internal_stops(protein))

writeXStringSet(AAStringSet(setNames(tr$protein, tr$id)), snakemake@output[["faa"]], width = 100000)
write_tsv(tr, snakemake@output[["tsv"]])

# For alignment: stops become X so mafft/iqtree accept them; references keep their own names (made safe)
ref_aln <- setNames(as.character(refs), safe_name(names(refs)))
rep_aln <- setNames(gsub("\\*", "X", tr$protein), tr$id)
stopifnot(!anyDuplicated(c(names(ref_aln), names(rep_aln))))
writeXStringSet(AAStringSet(c(ref_aln, rep_aln)), snakemake@output[["for_aln"]], width = 100000)
message(nrow(tr), " sequences translated; ", sum(tr$strand == "-"), " on the reverse strand; ",
        sum(tr$n_internal_stops > 0), " with internal stops")
