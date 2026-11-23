# Synthetic test dataset

A small, fully synthetic dataset to check an installation and to try the pipeline end to end in a few minutes.
No real samples are involved.

| File | Content |
|---|---|
| `source_sequences.fasta` | the six inserts (amplicon without primers) the reads are simulated from |
| `data/*.fastq.gz` | 3 samples x 400 simulated 2 x 300 bp MiSeq read pairs (about 420 KB in total) |
| `samples.tsv` | sample sheet with paths relative to this folder and a `Site` metadata column |
| `config.yaml` | test configuration (every read kept, the small test library below) |
| `subtype_library.tsv` | a 5-sequence curated library, so the test does not depend on `resources/subtype_library.tsv` |
| `subtype_assignments.tsv` | the true subtypes, standing in for the expert's assignments in the second run |

## What is in the samples

| Insert | Truth | TEST_01 | TEST_02 | TEST-03 | In the test library |
|---|---|---|---|---|---|
| `korv_A` | A | 200 | | 80 | yes |
| `korv_B` | B | 120 | | | yes |
| `korv_D` | D | 80 | 160 | | yes |
| `korv_L` | L | | 120 | 160 | yes |
| `korv_I_variant` | I | | 40 | 160 | no: 8 synonymous changes from a library I sequence (about 98.4% identity) |
| `nonfunctional_stop` | Non-functional | | 80 | | no: has an internal stop codon |

`TEST-03` has a dash on purpose, to show that sample names are made file-safe internally but reported unchanged.

Expected suggestions (one representative per insert, six in total):

| Insert | Suggested | Confidence | Flags |
|---|---|---|---|
| `korv_A`, `korv_B`, `korv_D`, `korv_L` | A, B, D, L | `high` | |
| `korv_I_variant` | I | `check` (98.4% to the library) | |
| `nonfunctional_stop` | Non-functional | `check` | `internal_stop`, `deletion_conserved_region` |

The read counts in the results are lower than in the table above because reads with sequencing errors form 100% OTUs with fewer than 3 reads, which are dropped, as in real data.
With the environments already built, both runs take about a minute.

## Run it

From the repository root, with the `korv` environment active:

```bash
snakemake --cores 2 --configfile .test/config.yaml
```

Open `results/korv_test/korv_test_subtype_review.html`, then run the report stage with the true assignments:

```bash
snakemake --cores 2 --configfile .test/config.yaml --config subtype_assignments=.test/subtype_assignments.tsv
```

This writes `results/korv_test/korv_test_final_result.xlsx`.

## Regenerating the reads

`tools/simulate_test_data.py` writes `data/` from `source_sequences.fasta` and `resources/ENV_primers.fasta`.
It is deterministic (fixed seed), so running it again gives identical files.
