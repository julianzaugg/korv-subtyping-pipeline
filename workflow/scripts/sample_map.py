# Sample sheet with the internal ids added, so reports show the original sample_id plus metadata.
import sys

sys.stderr = open(snakemake.log[0], "w")
rows = snakemake.params.rows
cols = list(rows[0].keys())
with open(snakemake.output[0], "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(r[c] for c in cols) + "\n")
