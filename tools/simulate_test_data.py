#!/usr/bin/env python3
"""Simulate the small paired-end MiSeq test dataset in .test/data/.

Each read pair comes from an amplicon  forward primer + insert + reverse-complemented
reverse primer, with the inserts taken from a FASTA file (.test/source_sequences.fasta)
and the primers from resources/ENV_primers.fasta. Reads are 2 x 300 bp with Phred
qualities that fall towards the 3' end and substitution errors drawn from them.
The output is deterministic for a given --seed.

The mix of inserts per sample is given as  sample_id:name=fraction,name=fraction ...

Usage (from the repository root, regenerates the committed test data):
  python3 tools/simulate_test_data.py
"""
import argparse
import gzip
import math
import random
from pathlib import Path

DEFAULT_MIX = [
    "TEST_01:korv_A=0.5,korv_B=0.3,korv_D=0.2",
    "TEST_02:korv_D=0.4,korv_L=0.3,nonfunctional_stop=0.2,korv_I_variant=0.1",
    "TEST_03:korv_A=0.2,korv_L=0.4,korv_I_variant=0.4",
]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--inserts", default=".test/source_sequences.fasta")
ap.add_argument("--primers", default="resources/ENV_primers.fasta")
ap.add_argument("--outdir", default=".test/data")
ap.add_argument("--pairs", type=int, default=400, help="read pairs per sample")
ap.add_argument("--read-length", type=int, default=300)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--mix", nargs="+", default=DEFAULT_MIX)
a = ap.parse_args()


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line and name:
            seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_from(template, rng):
    """One read with qualities falling from about Q36 to Q20 and matching substitution errors."""
    seq, qual = [], []
    n = len(template)
    for i, base in enumerate(template):
        q = max(2, min(40, round(36 - 16 * (i / n) ** 2 + rng.gauss(0, 2))))
        if rng.random() < 10 ** (-q / 10):
            base = rng.choice([b for b in "ACGT" if b != base])
        seq.append(base)
        qual.append(chr(q + 33))
    return "".join(seq), "".join(qual)


inserts = read_fasta(a.inserts)
fwd, rev = list(read_fasta(a.primers).values())
rng = random.Random(a.seed)
out = Path(a.outdir)
out.mkdir(parents=True, exist_ok=True)
L = a.read_length

for spec in a.mix:
    sample, parts = spec.split(":")
    mix = [(name, float(frac)) for name, frac in (p.split("=") for p in parts.split(","))]
    counts = {name: math.floor(frac * a.pairs) for name, frac in mix}
    counts[mix[0][0]] += a.pairs - sum(counts.values())
    order = [name for name, n in counts.items() for _ in range(n)]
    rng.shuffle(order)
    # Illumina-style names, so tools/make_sample_sheet.py also works on this folder
    paths = [out / f"{sample}_S1_L001_R{r}_001.fastq.gz" for r in (1, 2)]
    with gzip.GzipFile(paths[0], "wb", mtime=0) as g1, gzip.GzipFile(paths[1], "wb", mtime=0) as g2:
        for i, name in enumerate(order, 1):
            amplicon = fwd + inserts[name] + revcomp(rev)
            header = f"@SIM:{sample}:{i}:{name}"
            for g, mate, template in ((g1, 1, amplicon[:L]), (g2, 2, revcomp(amplicon)[:L])):
                s, q = read_from(template, rng)
                g.write(f"{header} {mate}:N:0:1\n{s}\n+\n{q}\n".encode())
    print(sample, ", ".join(f"{k}={v}" for k, v in counts.items()))
