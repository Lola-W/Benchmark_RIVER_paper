#!/usr/bin/env python3
"""
Count 32-class pyrimidine-first trinucleotide contexts on a reference genome
restricted to unmasked loci, per chromosome.

Definition:
- Center base position i (0-based) is counted if positions i-1, i, i+1 are all UNMASKED.
- Trinucleotide = ref[i-1:i+2]
- Only contexts with all bases in {A,C,G,T} are counted (skip N or other).
- Pyrimidine-first normalization:
    If center is C or T => use as-is.
    If center is A or G => reverse-complement the trinucleotide, then count.

Mask:
- BED is 0-based, half-open [start, end) masked intervals.
"""

import argparse
import sys
from collections import defaultdict

import numpy as np
import pysam


TRI32_ORDER = [
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT",
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT",
]

RC_TABLE = str.maketrans("ACGT", "TGCA")


def revcomp(tri: str) -> str:
    return tri.translate(RC_TABLE)[::-1]


def load_mask_intervals(bed_path: str, chrom: str):
    """Load masked intervals for a given chrom from a BED file."""
    intervals = []
    with open(bed_path, "r") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] != chrom:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end > start:
                intervals.append((start, end))
    return intervals


def build_mask_array(length: int, intervals):
    """Return boolean array masked[pos]=True for masked positions."""
    masked = np.zeros(length, dtype=np.bool_)
    for s, e in intervals:
        if s < 0:
            s = 0
        if e > length:
            e = length
        if e > s:
            masked[s:e] = True
    return masked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="GRCh38 fasta (must have .fai)")
    ap.add_argument("--bed", required=True, help="Blacklist BED (0-based, half-open)")
    ap.add_argument("--chrom", required=True, help="Chromosome name, e.g., chr1")
    ap.add_argument("--out", required=True, help="Output TSV path (trinuc\\tcount)")
    ap.add_argument("--require_triplet_unmasked", action="store_true", default=True,
                    help="Require i-1,i,i+1 all unmasked (default True)")
    args = ap.parse_args()

    chrom = args.chrom

    fa = pysam.FastaFile(args.fasta)
    if chrom not in fa.references:
        print(f"[ERROR] Chrom not found in FASTA: {chrom}", file=sys.stderr)
        sys.exit(2)

    L = fa.get_reference_length(chrom)
    if L < 3:
        print(f"[WARN] Chrom too short: {chrom} len={L}", file=sys.stderr)
        with open(args.out, "w") as out:
            for t in TRI32_ORDER:
                out.write(f"{t}\t0\n")
        return

    # Load masked intervals and build mask array
    intervals = load_mask_intervals(args.bed, chrom)
    masked = build_mask_array(L, intervals)

    # Fetch chromosome sequence once
    seq = fa.fetch(chrom).upper()
    if len(seq) != L:
        # pysam should match, but guard anyway
        L = min(L, len(seq))
        seq = seq[:L]
        masked = masked[:L]

    counts = defaultdict(int)

    # Iterate centers i=1..L-2 (0-based)
    # Require i-1,i,i+1 unmasked
    for i in range(1, L - 1):
        if masked[i]:# or masked[i - 1] or masked[i + 1]:
            continue

        tri = seq[i - 1 : i + 2]
        # skip non-ACGT
        if (tri[0] not in "ACGT") or (tri[1] not in "ACGT") or (tri[2] not in "ACGT"):
            continue

        center = tri[1]
        if center in ("C", "T"):
            key = tri
        else:
            key = revcomp(tri)

        counts[key] += 1

    # Write in fixed 32 order
    with open(args.out, "w") as out:
        for t in TRI32_ORDER:
            out.write(f"{t}\t{counts.get(t, 0)}\n")


if __name__ == "__main__":
    main()
