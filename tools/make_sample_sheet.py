#!/usr/bin/env python3
"""Draft a sample sheet (sample_id, R1, R2) from a folder of paired FASTQ files.

sample_id is the file name with the --strip regex removed (default: the Illumina
_S<n>_L<lane>_R<1|2>_001.fastq.gz tail). Add metadata columns to the result by hand.

Usage: make_sample_sheet.py FASTQ_DIR [--strip REGEX] [--recursive] > samples.tsv
"""
import argparse
import re
import sys
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("fastq_dir")
ap.add_argument("--strip", default=r"_S\d+_L\d+_R[12]_\d+\.f(ast)?q\.gz$")
ap.add_argument("--keep-sample-number", action="store_true",
                help="keep _S<n> in sample_id (e.g. SAMPLE01_S179) instead of stripping it")
ap.add_argument("--recursive", action="store_true")
a = ap.parse_args()

pat = re.compile(a.strip)
files = sorted(Path(a.fastq_dir).rglob("*.fastq.gz") if a.recursive else Path(a.fastq_dir).glob("*.fastq.gz"))
pairs = {}
for f in files:
    m = re.search(r"_R([12])_\d+\.f(ast)?q\.gz$", f.name)
    if not m:
        print(f"skipping {f.name}: no _R1/_R2 tag", file=sys.stderr)
        continue
    sid = pat.sub("", f.name)
    if a.keep_sample_number:
        sid = re.sub(r"_L\d+_R[12]_\d+\.f(ast)?q\.gz$", "", f.name)
    pairs.setdefault(sid, {})[m.group(1)] = f.resolve()

print("sample_id\tR1\tR2")
bad = 0
for sid, d in sorted(pairs.items()):
    if set(d) != {"1", "2"}:
        print(f"{sid}: missing a mate, skipped", file=sys.stderr)
        bad += 1
        continue
    print(f"{sid}\t{d['1']}\t{d['2']}")
sys.exit(1 if bad else 0)
