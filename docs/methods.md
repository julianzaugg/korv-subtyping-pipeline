# KoRV classification pipeline: methods

This document describes what the pipeline does and why, independent of any particular sequencing run.
Parameter values are the defaults in `config/config.yaml` and the files in `resources/`.
Software versions are those pinned in `workflow/envs/`.

## 1. Overview

The pipeline takes paired-end Illumina reads from an ENV amplicon of koala retrovirus (KoRV) and produces a per-sample table of KoRV subtype read counts and relative abundances.

It runs in four stages, orchestrated with Snakemake:

1. Preprocessing of each sample: merge read pairs, trim primers, filter, subsample, relabel.
2. Clustering of the pooled reads: 100% identity, abundance filter, then 97% identity.
3. Subtyping preparation: translation, alignment with reference subtypes, phylogeny, comparison with a curated library, and suggested subtypes with evidence.
4. Report: join the expert-assigned subtypes to the counts and write the final workbook.

The subtype of each sequence is always assigned by an expert.
The pipeline only makes suggestions and never writes a subtype into a result on its own.
Stage 4 does not run until the expert-completed assignments file is supplied.

## 2. Inputs

### 2.1 Sample sheet

A tab-separated file with one row per sample.
The required columns are `sample_id`, `R1` and `R2` (paths to the forward and reverse FASTQ files).
Any additional columns (for example a sample name, site or notes) are carried unchanged into the final report as sample metadata.
`sample_id` is the only key the pipeline uses.
It assumes no run number, LIMS identifier or file naming convention.

For use in read headers and file names, each `sample_id` is converted to an internal identifier by replacing every character other than letters, digits and underscores with `_`.
The pipeline stops if two samples share a `sample_id` or an internal identifier, or if a FASTQ file does not exist.
Reports always show the original `sample_id`.
`tools/make_sample_sheet.py` can draft a sheet from a folder of FASTQ files.

### 2.2 Reference resources

| File | Content |
|---|---|
| `resources/ENV_primers.fasta` | forward and reverse ENV primers, decoded from the original CLC adapter list |
| `resources/ENV_reference_sequence_translations_A-N.fasta` | 13 reference ENV protein sequences covering subtypes A to N |
| `resources/subtype_library.tsv` | curated expert assignments from 21 earlier jobs, one row per sequence per source (about 1,000 rows) |
| `resources/nonfunctional_rules.yaml` | thresholds for the non-functional flags and for the suggestion tiers |

The forward primer is `GCTTCTCATCTCAAACCCGC` and the reverse primer is `GGGTTGCCCGTAGGCGGT`.
The reference set is a fixed list of subtype representatives (for example KoRV-A, KoRV-C, KoRV-E, KoRV-L, KoRV-M, KoRV-N, D35, B21), not a database that grows with each run.
The library records past expert decisions for exact sequences and is plain text, so conflicting labels from different sources are kept and shown rather than resolved.

## 3. Stage 1: preprocessing (per sample)

All steps are run separately for each sample.

1. **Merging.**
   Read pairs are merged with `vsearch --fastq_mergepairs` (vsearch 2.28.1), with a minimum overlap of 20 bases (`merge_min_overlap`) and at most 6 mismatches in the overlap (`merge_max_diffs`).
   Pairs that do not merge are discarded.
2. **Primer trimming.**
   Both primers must be present, otherwise the read is discarded.
   This uses `cutadapt` 4.9 with a linked adapter: the forward primer at the 5' end and the reverse-complemented reverse primer at the 3' end, with `--discard-untrimmed`.
   Up to 12% error is allowed per primer (`primer_error_rate: 0.12`, which allows two mismatches in each primer).
3. **Quality, length and ambiguity filtering.**
   After the primers are removed, bases with quality below 13 are trimmed from both ends (`quality_cutoff`).
   Reads shorter than 400 bp or longer than 600 bp are discarded, as is any read containing an ambiguous base (`max_n: 0`).
   Primers are removed first so that end-quality trimming cannot cut into a primer and cause a good read to be lost.
4. **Subsampling.**
   If a sample has more than 15,000 reads after filtering, exactly 15,000 are drawn at random with `seqkit sample` using a fixed seed (`subsample_depth`, `subsample_seed`).
   Samples with 15,000 reads or fewer keep all of them.
   Subsampling can be switched off by setting the depth to `null`.
5. **Relabelling.**
   Reads are converted to upper-case FASTA and renamed `<internal_id>-<seven-digit read number>`, so the sample of origin travels in the read header.

A per-sample read statistics table records raw read pairs, merged reads and their average length, reads after trimming and filtering with the percentage retained and the average length, and reads after subsampling.

## 4. Stage 2: clustering

The reads of all samples are pooled and clustered together, so that the same sequence in two samples receives the same identifier.
The sample and read count of each sequence are tracked through vsearch `.uc` membership files and the count tables are rebuilt from those, instead of passing dense tables to the clustering tool.

1. **Dereplication.**
   Identical reads are collapsed with `vsearch --derep_fulllength` and their copy number is recorded.
2. **100% clustering.**
   The unique sequences are clustered at 100% identity with `vsearch --cluster_size`, using the copy numbers as abundance and with masking switched off.
   Because clustering at 100% identity also merges a sequence that is a prefix or substring of a longer one, this is not the same as simple dereplication.
   An OTU table of read counts per 100% OTU and sample is rebuilt from the dereplication and clustering `.uc` files.
3. **Abundance filter.**
   100% OTUs with fewer than 3 reads in total across all samples are removed (`min_total_reads_100`).
   This removes most sequencing-error variants.
   The number of reads removed per sample is recorded.
4. **Chimera removal.**
   The remaining 100% OTUs are checked with `vsearch --uchime3_denovo`, using their total read counts as abundance.
   OTUs called chimeric are removed from the tables and sequences before 97% clustering (`remove_chimeras`, default `true`), and the number of chimeric reads per sample is recorded.
   With `remove_chimeras: false` nothing is removed and the calls only appear as the `possible_chimera` flag.
5. **97% clustering.**
   The remaining 100% OTU sequences are clustered at 97% identity (`cluster_id_97`) with `vsearch --cluster_size`.
   In the default `legacy` order every 100% OTU is given abundance 1 when clustering, so each unique sequence counts equally, as in QIIME 2 style workflows.
   The alternative `abundance` order uses the real read counts, so that cluster centres are the dominant variants.
6. **Representative sequences.**
   The representative of each 97% cluster is its member 100% OTU with the most reads.
   Read counts of all members are summed per sample, giving a table of reads per representative and sample.
   The representatives are the sequences passed to Stage 3.

## 5. Stage 3: subtyping preparation

This stage produces evidence and suggestions for the expert reviewer.

### 5.1 Translation

Each representative is translated in the best of six reading frames (three on each strand).
"Best" is the frame whose protein gives the highest local alignment score (BLOSUM62, gap open 10, gap extend 1) against a subset of four reference proteins spread across the reference file.
No motif is required, so sequences that are damaged or divergent are still translated.
For the alignment, stop codons are written as `X`.

### 5.2 Alignment and tree

The representative proteins and the 13 references are aligned with MAFFT 7.526 (`--localpair --maxiterate 1000 --op 10 --ep 1`).
A maximum-likelihood tree is built with IQ-TREE 2.3.5, with ModelFinder restricted to the LG family of protein models (`-m MFP -mset LG`).

### 5.3 Chimera flag

The chimera check is run in Stage 2 (step 4).
By default chimeric OTUs are removed there, so none reach the representatives.
With `remove_chimeras: false` the representatives that were called chimeric are flagged `possible_chimera`, which is a flag only.

### 5.4 Comparison with the library

Each representative is searched against the curated library with `vsearch --usearch_global` (minimum identity 0.8, both strands, all candidates considered).
A representative whose oriented nucleotide sequence is identical to a library entry counts as an exact match.
The five closest library entries, their identity, labels and source job are reported.

### 5.5 Evidence per sequence

For each representative the pipeline reports:

- protein identity over aligned columns to every reference, the nearest reference and its subtype;
- the nearest reference of a different subtype and the identity margin to it (a small margin means the call is ambiguous);
- identity to the nearest reference over the conserved columns only (outside the NQQHNN zone and the variable region), which is what the novelty test uses;
- identity restricted to the variable region (the indel-rich block between the anchors `ASQQA` and `LGCSY` in the KoRV-A reference), which is reported for the expert (other than a region with fewer than 18 residues, indels there are not flagged);
- identity restricted to the NQQHNN zone near the start (between the anchors `GGAEM` and `PHQPM`), which varies slightly and tolerates some missing residues;
- the library evidence above;
- the number of internal stop codons;
- sequence missing at the 5' and 3' ends, deletions in conserved regions (excluding the NQQHNN zone and the variable region), missing sequence in the NQQHNN zone and in the variable region, in base pairs, measured against that sequence's own nearest reference, and the number of inserted residues.

### 5.6 Flags

A sequence is flagged when any of the following hold (thresholds in `resources/nonfunctional_rules.yaml`):

| Flag | Condition |
|---|---|
| `internal_stop` | any stop codon before the last residue |
| `deletion_conserved_region` | more than 4 bp deleted in a conserved region |
| `end_missing` | more than 4 bp missing at either end |
| `nqqhnn_deletion` | more than 5 bp missing in the NQQHNN zone |
| `variable_region_short` | fewer than 18 residues in the variable region (fewer than any complete reference) |
| `possible_chimera` | called by `uchime3_denovo` |

The first four are the criteria for suggesting "Non-functional".

### 5.7 Suggested subtype and confidence

The suggestion follows these rules in order:

1. If the best library hit has nucleotide identity of at least 97%, the library label is suggested.
   All hits within 1 percentage point of the best must carry the same label.
   The confidence is `high` if the identity is at least 99.5% and the labels agree.
   It is `check` if the identity is between 97% and 99.5%, or if the labels disagree (the majority label is then suggested).
   A `high` suggestion is lowered to `check` when non-functional features are present but the library label is not "Non-functional".
2. Otherwise, if any non-functional feature is present, the suggestion is "Non-functional" with confidence `check`.
   This comes before the novelty test because a damaged sequence is not a functional new subtype, however far it is from the references.
3. Otherwise, if the nearest library identity is below 95% (or there is no library hit) and the protein identity to the nearest reference over the conserved columns (outside the NQQHNN zone and the variable region) is below 80%, the suggestion is "possible new subtype (nearest X)" with confidence `novel`.
4. Otherwise, the subtype of the nearest reference is suggested with confidence `check`.

The basis of each suggestion is stated in words in the output.
`novel` means far from everything known, so these sequences always go to the expert.

### 5.8 Outputs for review

- `*_subtype_review.html`, a single self-contained page (no server or network needed) that shows the tree and the whole alignment with linked selection, the evidence for each sequence and controls to assign subtypes, and exports the completed workbook (see `docs/outputs.md`);
- `*_subtype_suggestions.xlsx`, with one row per representative, sorted with `novel` first and then `check` and `high`.
  It has blank `Subtype`, `Endogenous`, `Unsure` and `Note` columns for the expert, then the suggestion, confidence, basis, flags and all evidence columns.
  A second sheet explains each column.
- `review_figures/`:
  - two views of the tree of representatives and references: a cladogram (topology only) and a distance tree (branch lengths in expected amino acid substitutions per site), in which sequences with at least one flag are marked with a triangle;
  - a heatmap of protein identity to each reference, with the sequences ordered by hierarchical clustering of their identity profile and the references ordered by subtype;
  - an alignment of the variable region, cut from the full alignment, with rows in the order of the MAFFT alignment (guide-tree order, so similar sequences are neighbours) and cells coloured by difference from each sequence's nearest reference;
  - a single image of the whole alignment of representatives and references, with every sequence and column on one image however many sequences there are.
    Rows are in alignment order, references are in bold, the variable region is boxed and a top row shades each column by the share of its most common residue.
    Residues use the ClustalX colour scheme as described in the Jalview documentation: a residue is coloured only when its column meets that residue's composition criteria, and gaps, stop codons (`X`) and residues that do not meet the criteria stay white.
    The criteria percentages are of the residues in the column, gaps excluded.

  Each figure is written as PDF and PNG.
  The PDF is written with base `pdf()`, and the PNG with `ragg` when it is installed in the R environment.

## 6. Expert checkpoint

The expert reviews the suggestions, either in the interactive review page (tree and alignment side by side with linked selection) or in the suggestions workbook with the static figures.
The `Subtype` column is filled in for every representative (and `Endogenous` if used), either in the page, which exports the completed workbook, or directly in the workbook.
The file may be an .xlsx, .tsv or .csv.
The path is then given to the pipeline as `subtype_assignments` and the same command is run again.
Everything upstream is reused.
The pipeline stops with an error if any representative has no subtype.
Assignments are matched to representatives by `new_rep_97_ID`, or by sequence checksum (`seq_md5`) if the identifier is absent.

## 7. Stage 4: final report

The final workbook `*_final_result.xlsx` contains:

- read counts per sample and subtype, with the total number of KoRV reads per sample, and the relative abundance of each subtype (percentage of that sample's reads, rounded to two decimals), together with any sample metadata columns;
- the same table by subtype and endogenous status (`Subtype_Endogenous`), if an `Endogenous` column was filled in.
  The values are trimmed, spellings that differ only in case are merged, and blanks in a partly filled column are reported as `unassigned`.
  A column that is entirely blank is ignored;
- the number of samples in which each subtype occurs;
- the read counts before and after the abundance filter, and the per-sample read statistics from Stage 1;
- the expert's assignments with the representative sequences.

Samples with no KoRV reads after filtering have blank relative abundances.

## 8. Software

| Tool | Version | Used for |
|---|---|---|
| Snakemake | 9.27.2 (pinned in `environment.yml`, 8 or newer works) | workflow management and per-rule conda environments |
| vsearch | 2.28.1 | pair merging, dereplication, clustering, library search, chimera check |
| cutadapt | 4.9 | primer, quality and length filtering |
| seqkit | 2.8.2 | read statistics, subsampling, format conversion |
| MAFFT | 7.526 | protein alignment |
| IQ-TREE | 2.3.5 | maximum-likelihood tree and model selection |
| R | 4.5 | tables, translation, suggestions, figures and report (Biostrings, pwalign, ggtree, ggplot2, openxlsx and others) |

Primary references for the tools: vsearch (Rognes et al. 2016), cutadapt (Martin 2011), seqkit (Shen et al. 2016), MAFFT (Katoh and Standley 2013), IQ-TREE 2 (Minh et al. 2020), ModelFinder (Kalyaanamoorthy et al. 2017), UCHIME3 (Edgar 2016), ggtree (Yu et al. 2017).
Check the exact citations against the journals before publication.

## 9. Validation

The pipeline was validated against datasets that had already been analysed manually.

- **Clustering.**
  Starting from the reads of the manual workflow, the exact mode (no dereplication) reproduced the QIIME 2 results identically for two datasets.
  The default dereplicating mode differed by 0.07 to 0.14% of reads on three datasets, with the same 97% partitions in two of them.
- **Report.**
  From the manual workflow's reads, the report reproduced the curated workbook of one dataset exactly in exact mode, and within one read in three of twelve samples in the default mode.
- **From raw FASTQ.**
  On a 12-sample dataset, total KoRV reads were within 1% of the curated result and relative abundances within 0.25 percentage points.
  On a 60-sample dataset the median difference was 0.31 percentage points, with 51 of 60 samples within 2 points.
- **Suggestions.**
  With the dataset's own labels removed from the library, 78% of sequences (99.3% of reads) in one dataset received the expert's subtype, including 9 of 10 `high` suggestions.
  The misses were mainly sequences the expert had called non-functional without a flagged feature.

## 10. Limitations and caveats

- **Poor-quality runs.**
  The "at least 3 identical reads" filter is very sensitive to read noise.
  On a run with very poor reverse reads, vsearch's quality-aware merging gave cleaner reads than the manual workflow's merging, and so 18 to 55% more reads survived the filter and abundances differed by up to 10.8 percentage points from the manual result.
  Such samples will not reproduce manual numbers exactly.
- **Calibration.**
  The primer error rate and merge mismatch settings were tuned on two datasets and checked on a third.
- **Suggestion scoring.**
  Suggestions were scored on a small number of datasets, so the tier accuracy is a rough guide.
  The library is limited to what past experts assigned, and has some conflicting labels and mixed label vocabularies.
- **Scale.**
  100% clustering does not scale well.
  A 772,000-read dataset took 18 minutes, but a 432-sample dataset of about 1 GB of unique sequences took more than 1.5 hours.
  A faster variant that treats low-abundance sequences separately has been proposed but not implemented.
- **Subsampling.**
  Subsampling at 15,000 reads per sample is the default.
  Subsampling is random, so reruns with a different seed give slightly different low-abundance results.
  Choose `null` to keep all reads at the cost of more time and memory in clustering.
- **Platform.**
  The pipeline is tested on Linux and macOS (Apple Silicon); the synthetic test can be run on Linux on demand through GitHub Actions.
  The R steps use only the packages of their own conda environment and ignore any personal R library.

## 11. Reproducing a run

1. Install as described in the README, write the sample sheet and a config file (from `config/analysis_template.yaml`) giving `analysis_name`, `samples`, `outdir` and any changed parameters.
2. Run `snakemake --configfile <config> --cores <n> --conda-prefix <dir>`, which builds the conda environments on first use and writes `<outdir>/<analysis_name>/`, ending with the suggestions workbook and figures.
3. Review, fill in the assignments, add `subtype_assignments:` to the config and run the same command again to write the final workbook.

