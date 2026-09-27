#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys

VCF_HEADER = """##fileformat=VCFv4.2
##source=big_without_downsampling.annot.tsv
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

def parse_args():
    ap = argparse.ArgumentParser(description="TSV(chr,pos,ref,alt) -> N chunk VCFs (line-balanced). No GNU split needed.")
    ap.add_argument(
        "-i", "--input_tsv", required=True,
        help="Full path to big_without_downsampling.annot.tsv"
    )
    ap.add_argument(
        "-o", "--out_dir", required=True,
        help="Full path to output dir (will contain part_000.vcf ... part_099.vcf and chunk_list.txt)"
    )
    ap.add_argument(
        "-n", "--nchunks", type=int, default=100,
        help="Number of chunks (default: 100)"
    )
    return ap.parse_args()

def main():
    a = parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    # read TSV
    with open(a.input_tsv, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {k:i for i,k in enumerate(header)}
        for k in ["chr","pos","ref","alt"]:
            if k not in idx:
                sys.stderr.write(f"ERROR: missing column '{k}' in TSV header\n")
                sys.exit(2)

        rows = []
        for line in f:
            if not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            chrom = t[idx["chr"]]
            pos_s = t[idx["pos"]]
            ref = t[idx["ref"]]
            alt = t[idx["alt"]]
            if chrom == "." or pos_s == "." or ref == "." or alt == ".":
                continue
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            vid = f"{chrom}-{pos}-{ref}-{alt}"
            rows.append((chrom, pos, vid, ref, alt))

    if not rows:
        sys.stderr.write("ERROR: no valid variants parsed\n")
        sys.exit(2)

    # sort VCF-like
    rows.sort(key=lambda x: (x[0], x[1], x[3], x[4]))

    total = len(rows)
    nchunks = a.nchunks
    if nchunks < 1:
        sys.stderr.write("ERROR: nchunks must be >= 1\n")
        sys.exit(2)
    if nchunks > total:
        nchunks = total

    # line-balanced distribution
    base = total // nchunks
    rem = total % nchunks
    sizes = [base + (1 if i < rem else 0) for i in range(nchunks)]

    chunk_list = os.path.join(a.out_dir, "chunk_list.txt")
    with open(chunk_list, "w") as lf:
        idx_row = 0
        for i, sz in enumerate(sizes):
            out_vcf = os.path.join(a.out_dir, f"part_{i:03d}.vcf")
            with open(out_vcf, "w") as w:
                w.write(VCF_HEADER)
                for _ in range(sz):
                    chrom, pos, vid, ref, alt = rows[idx_row]
                    w.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t.\n")
                    idx_row += 1
            lf.write(out_vcf + "\n")

    sys.stderr.write(f"[OK] input={a.input_tsv}\n")
    sys.stderr.write(f"[OK] variants={total}\n")
    sys.stderr.write(f"[OK] out_dir={a.out_dir}\n")
    sys.stderr.write(f"[OK] nchunks={nchunks} (part_000.vcf .. part_{nchunks-1:03d}.vcf)\n")
    sys.stderr.write(f"[OK] chunk_list={chunk_list}\n")

if __name__ == "__main__":
    main()

