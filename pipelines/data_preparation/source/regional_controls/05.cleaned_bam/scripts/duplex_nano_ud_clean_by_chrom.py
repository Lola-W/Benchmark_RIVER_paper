#!/usr/bin/env python3
import argparse
import os
import pysam


def passes_read_qc(rec) -> bool:
# https://github.com/AlexandrovLab/DupCaller/blob/main/src/subcommands/funcs/call.py
    if rec.is_secondary:
        return False
    if rec.is_supplementary:
        return False

    if rec.is_qcfail:
        return False

    if not rec.is_proper_pair:
        return False

    if rec.has_tag("DT"):
        return False
    
    # call.py odd CIGAR first-op length==4 filter
    ct = rec.cigartuples
    if ct and len(ct) > 0:
        op, length = ct[0]
        if length == 4:
            return False

    return True


def canonical_umi_pair(rec) -> str:
    """
    Assumes qname ends with "..._UMI1+UMI2".
    Normalizes UMI order using (read1 forward) or (read2 reverse) => UMI1+UMI2 else swap.
    """
    bc = rec.query_name.rsplit("_", 1)[-1]
    if "+" not in bc:
        return bc
    bc1, bc2 = bc.split("+", 1)

    if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
        return f"{bc1}+{bc2}"
    else:
        return f"{bc2}+{bc1}"


def classify_f1r2_f2r1(rec) -> int:
    """
    Returns bit: F1R2=1, F2R1=2.
    F1R2: (read1 forward) OR (read2 reverse)
    F2R1: (read1 reverse) OR (read2 forward)
    """
    if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
        return 1
    return 2


def parse_region(region: str):
    # region like "chr1:1999800-2000500" (1-based inclusive for start)
    chrom, rest = region.split(":", 1)
    start_s, end_s = rest.split("-", 1)
    start0 = int(start_s.replace(",", "")) - 1
    end0 = int(end_s.replace(",", ""))
    return chrom, start0, end0


def main():
    ap = argparse.ArgumentParser(
        description="Create duplex-only BAM for a given chromosome (or region). Duplex key=(reference_start, canonical_UMI); duplex if >=1 F1R2 and >=1 F2R1."
    )
    ap.add_argument("--bam", required=True, help="Input BAM (indexed).")
    ap.add_argument("--chrom", default="", help="Chromosome name (e.g., chr1). Ignored if --region is provided.")
    ap.add_argument("--region", default="", help="Optional region like chr1:1999800-2000500.")
    ap.add_argument("--out-bam", required=True, help="Output BAM (unsorted; sbatch will sort+index).")
    ap.add_argument("--threads", type=int, default=1, help="Threads for BAM IO.")
    ap.add_argument("--progress-every", type=int, default=5_000_000, help="Progress print frequency.")
    args = ap.parse_args()

    if args.region:
        chrom, start0, end0 = parse_region(args.region)
        fetch_args = (chrom, start0, end0)
    else:
        if not args.chrom:
            raise SystemExit("ERROR: provide --chrom or --region")
        chrom = args.chrom
        fetch_args = (chrom,)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_bam)), exist_ok=True)

    bam_in = pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.threads))
    bam_out = pysam.AlignmentFile(args.out_bam, "wb", header=bam_in.header, threads=max(1, args.threads))

    # PASS 1: build duplex mask per (reference_start, canonical_umi)
    duplex_mask = {}  # key -> 0..3
    n_total = 0
    n_qc = 0

    for rec in bam_in.fetch(*fetch_args):
        n_total += 1
        if n_total % args.progress_every == 0:
            print(f"[PASS1] {chrom} scanned {n_total:,} reads; QC-pass {n_qc:,}; keys {len(duplex_mask):,}", flush=True)

        if not passes_read_qc(rec):
            continue
        n_qc += 1

        key = (rec.reference_start, canonical_umi_pair(rec))
        bit = classify_f1r2_f2r1(rec)
        prev = duplex_mask.get(key, 0)
        if prev != 3:
            duplex_mask[key] = prev | bit

    duplex_keys = sum(1 for v in duplex_mask.values() if v == 3)
    print(f"[PASS1 DONE] {chrom} total_reads={n_total:,} qc_pass={n_qc:,} keys={len(duplex_mask):,} duplex_keys={duplex_keys:,}", flush=True)

    # PASS 2: write only duplex reads
    n_written = 0
    for rec in bam_in.fetch(*fetch_args):
        if not passes_read_qc(rec):
            continue
        key = (rec.reference_start, canonical_umi_pair(rec))
        if duplex_mask.get(key, 0) == 3:
            bam_out.write(rec)
            n_written += 1
            if n_written % args.progress_every == 0:
                print(f"[PASS2] {chrom} written {n_written:,} reads", flush=True)

    bam_out.close()
    bam_in.close()

    print(f"[DONE] {chrom} wrote {n_written:,} reads (UNSORTED) -> {args.out_bam}", flush=True)


if __name__ == "__main__":
    main()
