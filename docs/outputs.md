# Pipeline outputs

This document explains every file the pipeline writes and every column of the suggestions workbook.
For how the values are computed see [methods.md](methods.md).

All files are in `results/<analysis_name>/` (the `outdir` of the config).
`<name>` below stands for `analysis_name`.

## 1. Which files do I need?

| You want to | Open |
|---|---|
| Assign subtypes (the expert step) | `<name>_subtype_review.html` (or edit `<name>_subtype_suggestions.xlsx` directly), with `review_figures/` for printing |
| Get the final counts per sample | `<name>_final_result.xlsx` (written after you supply the assignments) |
| Check read losses per sample | `<name>_sample_read_stats.tsv` |
| Re-use the sequences elsewhere | `<name>_final_sequences_to_subtype.fasta` and `.faa` |

Everything else is a record of intermediate steps.

## 2. Stage 1 and 2: reads and clustering

| File | Content |
|---|---|
| `<name>_sample_read_stats.tsv` | one row per sample: `Raw read pairs`; `Number of reads` and `Avg.length` (merged reads and their mean length); `Number of reads after trim` and `Avg.length after trim` (reads with both primers found that pass the quality, length and ambiguous-base filters, after primer removal); `Percentage trimmed (%)`, which despite its name is the percentage of merged reads that were kept; `Number of reads subsampled` (the reads used downstream, capped by `subsample_depth`) |
| `<name>_otu_table_100_filtered_counts.csv` | per sample: the number of reads in the 100% OTUs before and after dropping OTUs with fewer than 3 reads and chimeric OTUs, the number of reads removed by the read filter, and the number of reads in chimeric OTUs (`Reads_chimeric`) |
| `<name>_otu_table_100_greater_than_2.csv` | the abundance table of the 100% OTUs that passed the filter (one row per OTU, one column per sample, read counts) |
| `<name>_rep_sequences_100_greater_than_2.fasta` | the sequences of those 100% OTUs |
| `<name>_otu_table_100_with_reps.xlsx` | sheet `OTU_table_100_filt_with_reps`: each 100% OTU with the 97% cluster it belongs to and its counts per sample. Sheet `OTU_table_100_agg_per_rep_97`: the same counts summed per 97% representative (`new_rep_97_ID`, `Cluster_ID`, `N_samples_present`, `Total_number_of_reads`, then one column per sample) |
| `<name>_final_sequences_to_subtype.fasta` | the representatives to subtype, one per 97% cluster (nucleotide, in the orientation of the reads) |
| `<name>_final_sequences_to_subtype.tsv` | the same as a table: `new_rep_97_ID`, `Cluster_ID`, `N_samples_present`, `Total_number_of_reads`, `sequence`, `seq_md5` (checksum of the sequence, usable as a key) |
| `<name>_final_sequences_to_subtype.faa` | the representatives translated in their best frame |

The representative of a 97% cluster is the 100% member with the most reads.
Read headers are `<sample_id>-<number>`, so a representative ID such as `S1-0000010` only says which sample its first read came from.

## 3. Stage 3: subtyping preparation

| File | Content |
|---|---|
| `<name>_final_sequences_to_subtype_with_references.aln` | MAFFT protein alignment of the representatives and the A-N references (the row order is the MAFFT guide-tree order, so similar sequences are neighbours) |
| `<name>_final_sequences_to_subtype_with_references.treefile` | IQ-TREE maximum-likelihood tree of that alignment (also `.iqtree` report, `.log`, `.mldist`, `.bionj`, `.ckp.gz`, `.model.gz`) |
| `<name>_subtype_review.html` | the interactive review page (section 3.1) |
| `<name>_subtype_suggestions.xlsx` | suggestions, flags and evidence, with blank columns for you (section 4) |
| `review_figures/` | PDF and PNG of the figures in section 3.1 |
| `subtyping/translation.tsv` | the chosen translation of each representative: strand, frame, protein, alignment scores of the best and the runner-up frame, internal stop count |
| `subtyping/library_hits.tsv` | every library sequence at 80% nucleotide identity or more to each representative: query, library checksum, identity, query coverage, target coverage (no header) |
| `subtyping/chimeras.fasta` | the 100% OTUs called chimeric by `uchime3_denovo`; removed from the results when `remove_chimeras` is true (the default) |
| `subtyping/evidence.tsv` | the full evidence table behind the workbook, with extra columns (section 5) |
| `subtyping/logs/` | the logs of MAFFT, IQ-TREE, vsearch library search and uchime |

Also in `subtyping/`: `proteins_for_alignment.faa`, `final_sized.fasta` (sequences with read counts as abundance, for uchime) and `library.fasta` (the curated library as a search database).
These are inputs of the steps above.

### 3.1 The interactive review page

`<name>_subtype_review.html` is one self-contained file: open it in a browser by double clicking it, no server or internet connection is needed, and it can be emailed to a collaborator.
It shows the tree on the left and the whole alignment (ClustalX colours, variable region boxed, conservation row) on the right, with the evidence for the selected sequence and the assignment controls in a side panel.

- **Linked selection.**
  Click a tip, a branch (selects the whole clade), a row label or any alignment cell to select sequences, and they are highlighted in both views.
  Shift+click selects a range, Ctrl or Cmd+click toggles, and dragging over the labels selects a range.
- **Assigning.**
  Type a Subtype (plus the optional Endogenous, Unsure and Note) and press Enter, use a quick-label button, or accept the suggestion for the selected sequences.
  The quick labels never change order: keys 1 to 9 are the first nine subtypes of the references (A to I), 0 is Non-functional, and the remaining subtypes and any label you have used have unnumbered buttons.
  Sequences suggested as `possible new subtype` are never accepted in bulk: they need an explicit label.
  Labels that differ only in case, spaces, `-` or `_` are merged into the existing spelling, and any spelling of non-functional becomes `Non-functional`.
  You can also edit in the table: double-click the Assigned cell of a row, or press Enter to edit the selected rows, type a label and press Enter (Tab saves and moves to the next unassigned row, Esc cancels, an empty label clears it).
  Delete or Backspace clears the assignment of the selected rows.
  Undo and redo are available, and the top bar shows how many sequences still have no Subtype.
- **Looking around.**
  Rows can be ordered by suggested subtype (the default), assigned subtype, tree, alignment file, reads, workbook order or unassigned first.
  In the two subtype orders the rows are grouped by subtype, with the references of each subtype first and then the representatives in tree order, and a line marks each change of group.
  The groups follow the `subtype_order` section of `resources/nonfunctional_rules.yaml`: `order` lists the subtypes (alphabetical by default), and an intermediate such as QL or QD is placed next to its `nearer` parent, on the side facing its other parent, so it sits between the two as closely as one row order allows.
  A label written `X/Y` (such as A/D) is placed in the same way, next to X on the side facing Y.
  After the listed subtypes come `possible new subtype`, any other label, Non-functional and unlabelled rows.
  The assigned order uses the suggestion for sequences that have no Subtype yet.
  Rows are ordered when the order is chosen, so they do not jump while you assign.
  The tree can be a cladogram or a distance tree, tips can be coloured by the suggested or the assigned subtype, and the alignment can be coloured by ClustalX, by difference from the nearest reference, or not at all.
  The filters, which combine, are: `Show` (confidence tier), `has flags` (one or more flags, marked with a triangle), `unassigned` (no Subtype yet), `differs from suggestion` (the Subtype you assigned is not the suggested one, for reviewing your overrides), text search and minimum reads.
  They hide the rows that do not match (references stay visible), or dim them if "hide others" is unticked.
  Except in the tree order, or when rows are hidden, the tree keeps its own order and scroll position and the selection is still linked.
  A triangle marks a sequence with flags and a red ring marks a novel suggestion.
- **Saving.**
  Assignments are kept in the browser while you work, but browser storage can be cleared, so use **Export xlsx** to save.
  It writes `<name>_subtype_assignments.xlsx`: the suggestions workbook with the four columns filled in, which you pass to the pipeline as `subtype_assignments`.
  If any sequence has no Subtype the file is named `..._INCOMPLETE.xlsx` and the pipeline will refuse it until it is complete.
  **TSV** writes the same table as text, and **Session** saves a file that restores the assignments in the page.
- **Importing.**
  **Import** reads an exported workbook (including one you have edited in Excel), a TSV or CSV with a `new_rep_97_ID` (or `seq_md5`) column and a `Subtype` column, or a session file.
  Rows are matched by `new_rep_97_ID`, rows that are not in this run are reported and ignored, and where your assignments and the file disagree you choose which to keep.
  This is how to continue work in a later sitting, or to merge the calls of two reviewers.

The page needs a current browser (it uses canvas, and the browser's built-in decompression to read workbooks).
It was tested with up to about 400 sequences, where redrawing takes a few milliseconds.

### 3.2 Review figures

All figures have a PDF and a PNG version.
Rows and tips are labelled `<ID> [<suggested subtype>, <reads> reads, <n> samples]`.
References have a `REF` prefix (in the alignment figures) and are drawn in bold or as anchors.

| File | What it shows |
|---|---|
| `<name>_tree_cladogram` | tree topology only, so labels are easy to read |
| `<name>_tree_distance` | the same tree with branch lengths (expected amino acid substitutions per site) and a scale bar |
| `<name>_reference_identity_heatmap` | protein identity of each representative to each reference. Representatives are ordered by clustering their identity profile, references are grouped by subtype |
| `<name>_variable_region_alignment` | the variable region (between `ASQQA` and `LGCSY` in the KoRV-A reference) cut out of the full alignment, with four flanking columns on each side. Rows follow the alignment order. Orange cells differ from the sequence's nearest reference, white cells are gaps, blue rows are references |
| `<name>_alignment` | the whole alignment on a single image in ClustalX colours, with a top row showing how conserved each column is and the variable region boxed. Intended for pasting into the workbook you email |

In the tree figures a triangle marks a sequence with one or more flags (see section 4.3).

## 4. The suggestions workbook

`<name>_subtype_suggestions.xlsx` has two sheets: `suggestions` (one row per representative) and `legend` (a short column guide).
The first three columns are frozen.
Rows are sorted with `novel` first, then `check`, then `high`, then by suggested subtype and read count.
The three confidence labels are coloured (novel red, check orange, high green).
The pipeline only suggests: the call is yours.

### 4.1 What you fill in (yellow columns)

| Column | Content |
|---|---|
| `Subtype` | **Required for every row.** Your final subtype call. It is trimmed of spaces but otherwise used as typed, so spell each subtype the same way every time (`Non-functional`, not `non functional`) |
| `Endogenous` | Optional free text. Leave the whole column blank if you do not use it. If you fill it in, the final workbook gets an extra sheet by `Subtype_Endogenous`. Spaces are trimmed and spellings that differ only in case are merged (`yes` and `Yes`). Rows left blank in a partly filled column are reported as `unassigned`. The pipeline does not interpret the values (the curated library uses `Y`, `N` and `M`), so choose a vocabulary and keep to it |
| `Unsure` | Optional marker of confidence in the call, for you and your collaborators (for example `?`). It is kept in the final workbook's assignments sheet but is not used in any calculation |
| `Note` | Optional free text, kept in the assignments sheet only |

Do not change `new_rep_97_ID`: the assignments are matched to the sequences by it (or by `seq_md5` if the ID column is removed).

### 4.2 Identification

| Column | Content |
|---|---|
| `sample_id` | the sample with the most reads of this sequence. It is not a key and does not mean the sequence is absent elsewhere |
| `new_rep_97_ID` | identifier of the representative, `<sample>-<number>` |
| `seq_md5` | checksum of the nucleotide sequence, stable across runs and jobs |
| `Total_number_of_reads` | reads of this representative summed over all samples (all 100% members of its 97% cluster) |
| `N_samples_present` | number of samples where it has at least one read |

### 4.3 The suggestion

| Column | Content |
|---|---|
| `suggested_subtype` | the suggested label: a library label, `Non-functional`, a subtype letter from the nearest reference, or `possible new subtype (nearest X)` |
| `confidence` | `high`, `check` or `novel` (below) |
| `basis` | one line saying why, including the identity it rests on |
| `flags` | comma-separated list of criteria met (below), blank if none |

How the suggestion is made, in this order:

1. If a curated library sequence has at least 97% nucleotide identity (`check_library_identity`), the label is the one the best hits give.
   Hits within 1% identity of the best must agree: if they do, the tier is `high` at 99.5% identity or more and `check` between 97% and 99.5%.
   If they disagree, the majority label is used and the tier is `check`.
   A `high` is lowered to `check` if the sequence has non-functional features but the label is not non-functional.
2. Otherwise, if any non-functional feature is flagged: `Non-functional`, tier `check`.
   A damaged sequence is not a functional new subtype, so this comes before the next step even when the sequence is far from every reference (the basis then says so).
3. Otherwise, if the best library hit is below 95% identity (or there is none) and `conserved_identity` is below 80%: `possible new subtype`, tier `novel`.
   A new subtype is rare, so this is deliberately strict. In three previously assigned datasets none of the 123 functional sequences is below 85% in the conserved regions, while a 90% cut called 14 of them novel.
4. Otherwise the subtype of the nearest A-N reference, tier `check`.

Only a library match can give `high`.
All thresholds are in `resources/nonfunctional_rules.yaml`.

Flags:

| Flag | Raised when |
|---|---|
| `internal_stop` | a stop codon before the last residue |
| `deletion_conserved_region` | more than 4 bp missing inside the sequence, outside the NQQHNN zone and the variable region |
| `end_missing` | more than 4 bp missing at the 5' or the 3' end |
| `nqqhnn_deletion` | more than 5 bp missing in the NQQHNN zone |
| `variable_region_short` | fewer than 18 residues in the variable region, which is fewer than any complete reference (KoRV-A has 18). In the curated library 17 of the 18 such sequences without a stop codon are non-functional. In two previously assigned datasets it marks 43 of 60 non-functional sequences and none of the 102 with a functional label |
| `possible_chimera` | called by `uchime3_denovo` |

The first four and `variable_region_short` are the non-functional features that can lead to a `Non-functional` suggestion.
`possible_chimera` is shown but does not change the suggestion.
It only appears when `remove_chimeras` is false, because chimeric OTUs are otherwise removed before the representatives are chosen.
Other indels in the variable region are not flagged: they are expected and are for you to judge.

### 4.4 Nearest references

Identities are protein identities over the alignment columns where both sequences have a residue.
Missing sequence is counted in base pairs (alignment columns times three) against the sequence's own nearest reference.

| Column | Content |
|---|---|
| `nearest_reference` | the A-N reference with the highest protein identity over the whole alignment |
| `nearest_reference_identity` | that identity, from 0 to 1 |
| `conserved_identity` | the best identity to any reference over the conserved columns only, that is outside the NQQHNN zone and the variable region, where indels make identities lower for ordinary sequences. This is what the novelty test uses |
| `other_subtype_reference` | the best reference of a different subtype |
| `other_subtype_identity` | its identity |
| `margin_to_other_subtype` | `nearest_reference_identity` minus `other_subtype_identity`. A small margin means the call between two subtypes is ambiguous |

### 4.5 Variable region and NQQHNN zone

| Column | Content |
|---|---|
| `variable_region_nearest_reference` | the reference most similar within the variable region only (between `ASQQA` and `LGCSY`) |
| `variable_region_identity` | the identity to it within that region |
| `variable_region_missing_bp` | bp missing in the region. Informational only |
| `variable_region_residues` | number of residues the sequence has in the region (KoRV-A has 18). Below 18 raises `variable_region_short` |
| `nqqhnn_nearest_reference` | the reference most similar within the NQQHNN zone near the start (between `GGAEM` and `PHQPM`) |
| `nqqhnn_identity` | the identity to it within the zone |
| `missing_nqqhnn_bp` | bp missing in the zone. Flagged above 5 |

The variable region is where subtypes differ by indels, so the call there needs your judgment: compare it with the variable-region figure.

### 4.6 Library evidence

The library `resources/subtype_library.tsv` holds sequences with past expert assignments.
It is searched by nucleotide identity (both strands).

| Column | Content |
|---|---|
| `library_exact_match` | TRUE if this exact sequence is in the library |
| `library_best_identity` | nucleotide identity of the best hit, from 0 to 1, blank if no hit at 80% |
| `library_labels_best` | the label(s) the best hit carries. `A \| D` means two sources labelled that sequence differently |
| `library_top_hits` | up to five nearest library sequences as `label (identity %; source analysis)` |
| `library_endogenous` | the `Endogenous` value of the best hit, if it has one |

### 4.7 Sequence features

| Column | Content |
|---|---|
| `n_internal_stops` | stop codons before the last residue |
| `terminal_missing_5_bp`, `terminal_missing_3_bp` | bp missing at the start and at the end, relative to the nearest reference |
| `deletion_conserved_bp` | bp missing inside the sequence, outside the NQQHNN zone and the variable region |
| `inserted_residues` | residues the sequence has where its nearest reference has none |
| `strand`, `frame` | the strand (`+` or `-`) and reading frame (1 to 3) of the chosen translation |
| `Sequence` | the protein sequence, without the final stop |

## 5. Extra columns in `subtyping/evidence.tsv`

It has all the columns above except the four you fill in and `Sequence` (the sample column is called `top_sample` there), and in addition:

| Column | Content |
|---|---|
| `nearest_reference_subtype` | subtype letter of the nearest reference |
| `nearest_reference_coverage` | share of the nearest reference's residues that are aligned to a residue of the sequence |
| `library_labels_window`, `library_n_window_labels`, `library_majority_label` | the labels of the hits within 1% of the best, how many distinct labels there are, and the most common one (what the suggestion uses) |
| `chimera` | TRUE if called chimeric |
| `nonfunctional_features` | the flags that count as non-functional features |
| `protein` | the protein sequence with the stop marked `*` |
| `identity__<reference>` | the protein identity to every reference, one column each |

## 6. Stage 4: the final workbook

Written by the second run, once `subtype_assignments` is set.
The run stops with an error if any representative has no `Subtype`.

All sheets use Arial, centred cells, bold column names and column widths fitted to the content.
On the counts sheets the relative abundances are shaded from white (0%) to dark red (100%).

| Sheet | Content |
|---|---|
| `counts_and_abundances` | one row per sample: your sample sheet's extra columns (for example `Name`), `Total_KoRV_reads`, then read counts per subtype, then relative abundance per subtype (percent of that sample's KoRV reads, two decimals, blank if the sample has no reads). The last row counts the samples where each subtype is present |
| `counts_and_abundances_subtype`, `counts_and_abundances_sub_endo` | the same table by `Subtype`, and by `Subtype_Endogenous` (for example `A_N`). Written instead of the sheet above only if `Endogenous` was filled in |
| `Reads_filtered_OTU_100` | the content of `*_otu_table_100_filtered_counts.csv` |
| `sample_read_stats` | the content of `*_sample_read_stats.tsv` |
| `Representative_KorV_sequences` | your assignments sheet as supplied, so the call and its evidence stay together |

## 7. The `work/` directory

`work/` holds the intermediate files (per-sample merged and filtered reads, the vsearch clustering files and the read-by-OTU table) and is not needed to read any result.
It is kept so a later run can reuse the earlier steps, and can be deleted once the analysis is final.
