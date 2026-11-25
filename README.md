# KoRV classification pipeline

Raw paired MiSeq reads of the KoRV ENV amplicon to a final KoRV subtype abundance workbook, without CLC or QIIME.
The methods are described in [docs/methods.md](docs/methods.md) and every output in [docs/outputs.md](docs/outputs.md).

## What it does

1. **Preprocess** (per sample): merge read pairs, trim the ENV primers (both required), filter by length and ambiguous bases, optionally subsample, relabel reads `<sample_id>-<number>`.
2. **Cluster**: dereplicate, cluster at 100%, drop OTUs with fewer than 3 reads, cluster at 97%, choose the most abundant 100% member as each 97% representative.
3. **Subtyping preparation**: translate each representative in its best frame, align with the A-N references, build a tree, compare with a curated library of past expert assignments, and write suggestions with evidence, review figures and an interactive review page.
4. **Report**: join the expert-completed assignments and write the final workbook.

The subtype call itself is always made by the expert.
The pipeline only suggests, and the report stage needs the completed assignments.

## Installation

The pipeline runs on Linux, macOS (Intel and Apple Silicon) and Windows through WSL.
It needs only conda: Snakemake runs in one small environment, and every step installs its own tools (vsearch, cutadapt, seqkit, MAFFT, IQ-TREE, R and packages) the first time it runs.
No administrator rights are needed.

1. Install conda if you do not have it.
   [Miniforge](https://conda-forge.org/download/) is recommended; any conda 24.7.1 or newer, or mamba, also works.
   Set strict channel priority once, which conda-forge recommends and Snakemake warns about otherwise:

   ```bash
   conda config --set channel_priority strict
   ```

2. Get the code:

   ```bash
   git clone https://github.com/julianzaugg/korv-subtyping-pipeline.git
   ```

   ```bash
   cd korv_classification_pipeline
   ```

3. Create and activate the Snakemake environment ([environment.yml](environment.yml)):

   ```bash
   conda env create -f environment.yml
   ```

   ```bash
   conda activate korv
   ```

4. Run the test (next section).
   The first run builds the step environments, which takes about 5 minutes and 1.5 GB of disk (longer on a slow connection); later runs reuse them.

Run all `snakemake` commands from the repository folder.
By default the step environments are stored in `.snakemake/conda/` inside the clone.
To share them between clones or projects, add `--conda-prefix ~/.korv_conda_envs` (any folder) to every command.

## Quick test

A small synthetic dataset (3 samples, 1,200 read pairs, six known sequences) is included in [.test/](.test/README.md).

```bash
snakemake --cores 2 --configfile .test/config.yaml
```

This writes `results/korv_test/`.
Open `results/korv_test/korv_test_subtype_review.html` in a browser to see the review page.
Then run the report stage with the known answers standing in for the expert:

```bash
snakemake --cores 2 --configfile .test/config.yaml --config subtype_assignments=.test/subtype_assignments.tsv
```

If `results/korv_test/korv_test_final_result.xlsx` exists, the installation works.
[.test/README.md](.test/README.md) lists what the suggestions should be.

## Running your own data

### 1. Sample sheet

A tab-separated file with the required columns `sample_id`, `R1` and `R2`; see [config/samples_template.tsv](config/samples_template.tsv).
Any other columns (for example `Name`, `Site`) are carried through to the report as metadata.
Relative `R1`/`R2` paths are read from the folder that holds the sample sheet.
`sample_id` is the only key the pipeline uses.
Characters other than letters, digits and underscores are replaced by `_` in read headers and file names (and checked for uniqueness), but reports always show the original `sample_id`.

To draft one from a folder of Illumina-named FASTQ files (`<name>_S1_L001_R1_001.fastq.gz`):

```bash
python tools/make_sample_sheet.py /path/to/fastq_dir > /path/to/my_analysis/samples.tsv
```

### 2. Configuration

Copy [config/analysis_template.yaml](config/analysis_template.yaml) next to your data and set `analysis_name`, `samples` and `outdir`.
Every other setting has a default in [config/config.yaml](config/config.yaml) (filters, clustering, thresholds), which you can override by adding the same key to your file.
A misspelled key stops the run with an error.

### 3. Suggestions

```bash
snakemake --cores 8 --configfile /path/to/my_analysis/config.yaml
```

`--cores` is the total number of CPU cores Snakemake may use; `threads` in the config (default 4) is what a single multi-threaded step uses.
On a shared machine set `--cores` to what you are allowed to use.
Use `-n` (dry run) to see what would run without running it.

### 4. Expert assignment and report

Assign a subtype to every representative sequence in either of two ways:

- **Review page**: open `<outdir>/<analysis_name>/<analysis_name>_subtype_review.html` in a browser (no internet or server needed), assign the subtypes and click Export xlsx.
- **Workbook only**: open `<outdir>/<analysis_name>/<analysis_name>_subtype_suggestions.xlsx` in Excel (or LibreOffice) and fill in the `Subtype` column, plus `Endogenous`, `Unsure` and `Note` if you use them.
  Save it under a new name so a later run cannot overwrite your work.
  The figures in `review_figures/` show the same evidence as the page.

Both give a workbook in the same format, and you can mix them (the page imports an edited workbook).
Set `subtype_assignments: /path/to/the/completed.xlsx` in your config and run the same command again.
The second run reuses everything upstream and writes `<analysis_name>_final_result.xlsx`.
It refuses to finish if any representative has no subtype.

### Running on an HPC cluster

The simplest way needs nothing extra: submit the normal command as one batch job on one node, and give `--cores` the number of cores the job asked for.
Snakemake then runs the steps inside that job exactly as on a laptop.
One node is enough because the slow step (100% clustering) is a single job and cannot be split across nodes anyway.
A dataset of a few hundred thousand reads needs about 4 cores, 4 GB of memory and well under an hour.

1. Compute nodes often have no internet access, so build the environments once on a login node:

   ```bash
   snakemake --configfile /path/to/my_analysis/config.yaml --conda-create-envs-only --cores 1
   ```

2. Submit a job script like this one (SLURM shown; for PBS use `#PBS -l select=1:ncpus=8:mem=16gb` and `#PBS -l walltime=04:00:00`, and so on for other schedulers):

   ```bash
   #!/bin/bash
   #SBATCH --job-name=korv
   #SBATCH --nodes=1
   #SBATCH --cpus-per-task=8
   #SBATCH --mem=16G
   #SBATCH --time=04:00:00

   conda activate korv
   cd /path/to/korv_classification_pipeline
   snakemake --cores 8 --configfile /path/to/my_analysis/config.yaml
   ```

   If `conda activate` fails inside a job, add `eval "$(conda shell.bash hook)"` before it, or load your site's conda module.

Snakemake can also submit every step as a separate cluster job, through an [executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/) for your scheduler (for example `snakemake-executor-plugin-slurm`).
That is an extra install and needs site-specific settings (account, partition), and it only helps with many samples, because only the per-sample preprocessing runs in parallel.
For this pipeline the single job above is recommended.

## Main outputs (`results/<analysis_name>/`)

Every file and every column of the suggestions workbook is explained in [docs/outputs.md](docs/outputs.md).

| File | Content |
|---|---|
| `*_sample_read_stats.tsv` | reads per processing step per sample (same columns as the old `*_data.xlsx`) |
| `*_otu_table_100_filtered_counts.csv` | reads before and after the minimum-read and chimera filters per sample |
| `*_otu_table_100_with_reps.xlsx` | 100% OTUs with their 97% cluster, and counts collapsed per 97% representative |
| `*_final_sequences_to_subtype.fasta/.tsv/.faa` | the representatives to subtype (nucleotide, table with `seq_md5`, protein) |
| `*_final_sequences_to_subtype_with_references.aln/.treefile/.iqtree` | alignment and tree with the A-N references |
| `*_subtype_review.html` | interactive page: tree and whole alignment side by side with linked selection, evidence per sequence, and subtype assignment that exports the completed workbook (see [docs/outputs.md](docs/outputs.md#31-the-interactive-review-page)) |
| `*_subtype_suggestions.xlsx` | suggestions, confidence tier, flags, evidence and blank columns for the expert |
| `review_figures/` | cladogram and distance tree, reference identity heatmap, variable-region alignment, whole alignment in ClustalX colours (PDF and PNG) |
| `*_final_result.xlsx` | counts and relative abundances per sample, plus supporting sheets |

## Suggestions, in short

The workbook has a lot of columns.
You fill in the four yellow ones (`Subtype` is required) and the rest is evidence: see [docs/outputs.md](docs/outputs.md#4-the-suggestions-workbook).

- Translation uses the best of six frames (no motif is required), checked against the previously curated proteins of all 503 functional-subtype sequences in the library.
- `high`: a close match (at least 99.5% nucleotide identity) to curated library sequences that agree.
  `check`: needs a look (library labels disagree, lower identity, or flags).
  `novel`: rare, below 80% identity to every A-N reference in the conserved regions and no close library sequence.
- Flags follow `resources/nonfunctional_rules.yaml`: internal stop codons, deletions in conserved regions (over 4 bp), sequence missing at the ends, missing residues in the NQQHNN zone (over 5 bp), possible chimeras.
  Indels in the variable region (between `ASQQA` and `LGCSY`) are not flagged, they are shown in the evidence columns and the figures for the expert to judge. The one exception is a region with fewer than 18 residues (`variable_region_short`), which is almost always non-functional.
  All thresholds are editable.
- The library `resources/subtype_library.tsv` was built from 21 legacy jobs by `tools/build_subtype_library.R`.
  It is a plain TSV (one row per sequence per source) and conflicts between sources are shown, not resolved.

## Troubleshooting

- **`CreateCondaEnvironmentException` or a message that conda is too old**: Snakemake needs conda 24.7.1 or newer.
  The `korv` environment includes a recent conda; make sure it is activated (`which snakemake` should point into it).
- **A step fails**: the message names its log file (under `<outdir>/<analysis_name>/work/logs/` or `subtyping/logs/`).
  Re-running the same command resumes from the failed step.
- **`ValidationError ... 'false' is not of type 'boolean'`**: `--config` passes `true`/`false` as text.
  Set yes/no options such as `dereplicate` in the config file instead.
- **`LockException`**: an earlier run was interrupted. Run the same command once with `--unlock`, then again without it.
- **`No such file or directory` for `config/config.yaml`**: run `snakemake` from the repository folder.
- **Apple Silicon**: all tools have native builds and the test passes natively.
  If a package is ever missing for `osx-arm64`, create the environments for Intel and let Rosetta run them by prefixing the command with `CONDA_SUBDIR=osx-64`.
- **Personal R packages**: the R steps ignore any R library of your own (for example `~/Library/R/...` from a CRAN R install) and use only their conda environment, so an existing R setup cannot interfere.
- **Windows**: use WSL (Ubuntu) and keep data inside the WSL file system for speed.

## Layout

- `workflow/` Snakefile, rules, conda environments, R scripts, config schema
- `config/` default settings, analysis template, sample sheet template
- `resources/` A-N reference translations, ENV primers, curated subtype library, non-functional rules
- `tools/` stand-alone helpers (sample sheet drafting, test data simulation, library builder, CLC adapter extraction)
- `.test/` synthetic test dataset
- `tests/` unit tests for the review page (`review_page/`) and validation scripts against legacy jobs (`legacy/`)
- `docs/` methods and outputs

## Development and validation

- `snakemake --lint --configfile .test/config.yaml` checks the workflow code, and `deno test --allow-read --allow-env tests/review_page/core_test.js` tests the review page logic.
  Both, and the synthetic test on Linux, can be run on GitHub Actions on demand ([.github/workflows/test.yml](.github/workflows/test.yml)): Actions tab, then Run workflow, or `gh workflow run test.yml`.
- `tools/simulate_test_data.py` regenerates `.test/data/` (deterministic).
- `start_from: subsetted_fasta` starts from a legacy CLC `*_subsetted.fasta` (headers `<sample>-<number>`) and `dereplicate: false` reproduces the old QIIME behaviour exactly.
  The scripts in `tests/legacy/` compare outputs with legacy results by sequence content; they need the legacy job outputs, which are not in this repository.

## Provenance of the primers

`resources/clc/ENV adapters.clc` is the original CLC adapter list.
`tools/extract_clc_adapters.py` decodes it to `resources/ENV_primers.fasta`.
