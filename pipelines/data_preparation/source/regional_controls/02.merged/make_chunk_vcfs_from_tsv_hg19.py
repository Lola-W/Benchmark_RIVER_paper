#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import sys
from typing import Dict, List, Tuple

VCF_HEADER = """##fileformat=VCFv4.2
##source=tsv_with_hg19_columns
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

def parse_args():
    ap = argparse.ArgumentParser(
        description="TSV(with hg19 columns) -> N chunk VCFs, preserving ORIGINAL TSV row order and chunk boundaries. "
                    "Each VCF record uses hg19_chr/pos/ref/alt, and stores RID + hg38 coords in INFO."
    )
    ap.add_argument("-i", "--input_tsv", required=True, help="Input TSV containing hg19_* columns")
    ap.add_argument("-o", "--out_dir", required=True, help="Output dir (part_000.vcf .., chunk_list.txt)")
    ap.add_argument("-n", "--nchunks", type=int, default=100, help="Number of chunks (default: 100)")
    return ap.parse_args()

def get_col_index(header: List[str], name: str) -> int:
    if name not in header:
        sys.stderr.write(f"ERROR: missing column '{name}' in TSV header\n")
        sys.exit(2)
    return header.index(name)

def is_valid_allele(a: str) -> bool:
    # English comment: allow A/C/G/T/N only, and non-empty.
    if not a or a == ".":
        return False
    for c in a:
        if c not in "ACGTNacgtn":
            return False
    return True

def main():
    args = parse_args()
    if args.nchunks < 1:
        sys.stderr.write("ERROR: --nchunks must be >= 1\n")
        sys.exit(2)

    os.makedirs(args.out_dir, exist_ok=True)

    # -----------------------------
    # Read header + stream body
    # -----------------------------
    with open(args.input_tsv, "r") as f:
        header = f.readline().rstrip("\n").split("\t")

        # Required columns
        i_chr38 = get_col_index(header, "chr")
        i_pos38 = get_col_index(header, "pos")

        i_chr19 = get_col_index(header, "hg19_chr")
        i_pos19 = get_col_index(header, "hg19_pos")
        i_ref19 = get_col_index(header, "hg19_ref")
        i_alt19 = get_col_index(header, "hg19_alt")

        # Optional: hg19_unmapped
        i_unmap = header.index("hg19_unmapped") if "hg19_unmapped" in header else -1

        # First pass: count total body rows (for chunk boundaries)
        # English comment: we need total rows INCLUDING unmapped placeholders to preserve original boundaries.
        body_lines = f.readlines()

    total_rows = len(body_lines)  # includes unmapped/placeholder rows
    if total_rows == 0:
        sys.stderr.write("ERROR: no body rows in TSV\n")
        sys.exit(2)

    nchunks = min(args.nchunks, total_rows)
    base, rem = divmod(total_rows, nchunks)
    sizes = [base + (1 if i < rem else 0) for i in range(nchunks)]

    # Precompute chunk boundaries on RID (1-based)
    # chunk i covers RID in [start, end] inclusive
    bounds: List[Tuple[int, int]] = []
    rid = 1
    for sz in sizes:
        start = rid
        end = rid + sz - 1
        bounds.append((start, end))
        rid = end + 1

    # Create writers
    out_paths: List[str] = []
    writers = []
    try:
        for i in range(nchunks):
            out_vcf = os.path.join(args.out_dir, f"part_{i:03d}.vcf")
            out_paths.append(out_vcf)
            w = open(out_vcf, "w")
            w.write(VCF_HEADER)
            writers.append(w)

        # Write chunk_list.txt
        chunk_list = os.path.join(args.out_dir, "chunk_list.txt")
        with open(chunk_list, "w") as lf:
            for p in out_paths:
                lf.write(p + "\n")

        # Helper: find chunk index from RID using bounds (linear scan is fine for 100 chunks)
        def chunk_of_rid(r: int) -> int:
            for i, (s, e) in enumerate(bounds):
                if s <= r <= e:
                    return i
            return nchunks - 1  # fallback (should not happen)

        # Second pass: iterate original order, assign by RID, write mapped-only records
        for rid, line in enumerate(body_lines, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            t = line.split("\t")
            ci = chunk_of_rid(rid)
            w = writers[ci]

            # Determine unmapped
            unmapped = False
            if i_unmap != -1:
                unmapped = (t[i_unmap] == "1")
            else:
                unmapped = (t[i_chr19] == "." or t[i_pos19] == ".")

            if unmapped:
                continue  # keep placeholder in boundaries, but no VCF record

            chr19 = t[i_chr19]
            pos19_s = t[i_pos19]
            ref19 = t[i_ref19]
            alt19 = t[i_alt19]

            # Basic sanity
            if chr19 in (".", "") or pos19_s in (".", ""):
                continue
            try:
                pos19 = int(pos19_s)
            except ValueError:
                continue
            if pos19 < 1:
                continue

            if not is_valid_allele(ref19) or not is_valid_allele(alt19):
                continue

            chr38 = t[i_chr38]
            pos38_s = t[i_pos38]

            # INFO fields: RID + hg38 coords
            info = f"RID={rid};HG38CHR={chr38};HG38POS={pos38_s}"

            # VCF: CHROM POS ID REF ALT QUAL FILTER INFO
            # English comment: ID is RID for stable join.
            w.write(f"{chr19}\t{pos19}\t{rid}\t{ref19}\t{alt19}\t.\tPASS\t{info}\n")

    finally:
        for w in writers:
            try:
                w.close()
            except Exception:
                pass

    sys.stderr.write(f"[OK] input={args.input_tsv}\n")
    sys.stderr.write(f"[OK] total_rows(body)={total_rows}\n")
    sys.stderr.write(f"[OK] nchunks={nchunks} (part_000.vcf .. part_{nchunks-1:03d}.vcf)\n")
    sys.stderr.write(f"[OK] out_dir={args.out_dir}\n")

if __name__ == "__main__":
    main()

