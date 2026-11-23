# korv_classification_pipeline

Snakemake + R pipeline replacing the manual CLC / QIIME / R workflow for KoRV subtyping.
User docs: README.md, docs/methods.md, docs/outputs.md.

## Rules for this repo
- Pipeline code must not assume job numbers, LIMS IDs, naming conventions or a particular machine or cluster. A sample is `sample_id` plus two FASTQ files.
- No hard-coded paths or `setwd()`. R scripts take paths from Snakemake.
- Final subtype assignment is manual (expert). The pipeline only suggests.
- Raw reads, qza/biom files and xlsx data never go in git. Never commit credentials. The only data in git is the synthetic test set in `.test/`.
- After a change, run the synthetic test (`snakemake --cores 2 --configfile .test/config.yaml`, see .test/README.md); where legacy jobs are available, also validate against them (tests/legacy/).
- Do not use em dashes in text. Commit messages carry no co-author lines.
