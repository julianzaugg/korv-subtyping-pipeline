# Helper functions for the Snakefile and the rules (kept apart from the rules, as Snakemake recommends).


def internal_id(sample_id):
    """File- and header-safe form of a sample_id. Reports always show the original."""
    return re.sub(r"[^A-Za-z0-9_]", "_", sample_id)


def load_samples(path):
    """Read the sample sheet. Relative R1/R2 paths are relative to the sample sheet's folder."""
    base = Path(path).parent
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no samples")
    for col in ("sample_id", "R1", "R2"):
        if col not in rows[0]:
            raise ValueError(f"{path}: missing required column '{col}'")
    seen = {}
    for r in rows:
        r["internal_id"] = internal_id(r["sample_id"])
        if r["sample_id"] in seen or r["internal_id"] in seen.values():
            raise ValueError(f"{path}: duplicate sample_id or internal id for '{r['sample_id']}'")
        seen[r["sample_id"]] = r["internal_id"]
        for col in ("R1", "R2"):
            f = Path(r[col]).expanduser()
            r[col] = str(f if f.is_absolute() else base / f)
            if not Path(r[col]).exists():
                raise ValueError(f"{path}: {r['sample_id']}: {col} not found: {r[col]}")
    return rows


def read_fasta_seqs(path):
    seqs, cur = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = []
            seqs.append(cur)
        elif line and cur is not None:
            cur.append(line)
    return ["".join(x).upper() for x in seqs]


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def raw_reads(wildcards):
    r = BY_INTERNAL[wildcards.sid]
    return {"r1": r["R1"], "r2": r["R2"]}


def report_inputs(wildcards):
    files = {
        "assignments": config["subtype_assignments"],
        "collapsed": f"{OUT}/{NAME}_otu_table_100_with_reps.xlsx",
        "final": f"{OUT}/{NAME}_final_sequences_to_subtype.tsv",
        "filtered_counts": f"{OUT}/{NAME}_otu_table_100_filtered_counts.csv",
    }
    if START == "fastq":
        files["sample_map"] = f"{WORK}/sample_map.tsv"
        files["read_stats"] = f"{OUT}/{NAME}_sample_read_stats.tsv"
    return files
