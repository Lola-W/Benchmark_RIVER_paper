#!/usr/bin/env python3
"""
vcf_to_vcf_overlap_filter.py

Purpose:
    Remove variants from an input VCF that also appear in a second
    ("filter") VCF. Variants are matched by exact CHROM, POS, REF and ALT.

Usage:
    python vcf_to_vcf_overlap_filter.py \
        --input-vcf  <input.vcf[.gz]> \
        --filter-vcf <filter.vcf[.gz]> \
        --output-vcf-gz <output.vcf.gz>

Inputs:
    --input-vcf      VCF to be filtered (plain or gzipped)
    --filter-vcf     VCF listing variants to remove (plain or gzipped)

Output:
    --output-vcf-gz  Filtered VCF, bgzipped and tabix-indexed (.tbi)

Notes:
    - For multi-allelic records, the whole record is removed if ANY of its
      ALT alleles matches an entry in the filter VCF.
    - CHROM naming (e.g. "chr1" vs "1") must be identical in both files.
    - The filter VCF is held in memory as a set, so memory use scales with
      the number of filter variants (relevant for gnomAD).

Pipeline context:
    Snakemake rules step01_germline, step05_pon and step06_gnomad
    (filter VCF = germline calls, panel of normals and gnomAD, respectively).
"""

import argparse
import gzip
import subprocess
import sys
import tempfile
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter VCF by exact CHROM-POS-REF-ALT overlap against another VCF."
    )
    parser.add_argument("--input-vcf",    required=True, help="Input VCF to be filtered (.vcf or .vcf.gz).")          # Add input-vcf argument
    parser.add_argument("--filter-vcf",   required=True, help="VCF containing variants to remove (.vcf or .vcf.gz).") # Add filter-vcf argument
    parser.add_argument("--output-vcf-gz",required=True, help="Output filtered VCF.gz. Tabix index will also be created.") # Add output-vcf-gz argument
    return parser.parse_args()


def open_vcf(path):
    """Open .vcf or .vcf.gz as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")  # Open as gzip-compressed text file
    return open(path, "rt")           # Open as text file


def variant_keys(fields):
    """Return (CHROM, POS, REF, ALT) keys for all ALT alleles in a VCF record."""
    if len(fields) < 5:
        return []
    chrom, pos, ref = fields[0], fields[1], fields[3]
    return [
        (chrom, pos, ref, alt.strip())
        for alt in fields[4].split(",")
        if alt.strip() and alt.strip() != "."  # Add (CHROM, POS, REF, ALT) tuple as a key
    ]


def load_filter_keys(filter_vcf):
    """Load all CHROM-POS-REF-ALT keys from filter VCF into a set."""
    filter_keys = set()
    n_records = 0
    with open_vcf(filter_vcf) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            keys = variant_keys(line.rstrip("\n").split("\t"))
            if keys:
                n_records += 1
                filter_keys.update(keys)
    return filter_keys, n_records


def filter_vcf_to_temp(input_vcf, filter_keys, temp_vcf):
    """Filter input VCF against filter_keys, write to temporary uncompressed VCF."""
    n_total = n_removed = n_kept = 0
    with open_vcf(input_vcf) as fin, open(temp_vcf, "wt") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            n_total += 1
            # Remove the record if any of its ALT alleles is in the filter set
            if any(key in filter_keys for key in variant_keys(fields)):
                n_removed += 1
                continue  # Skip writing to output
            fout.write(line + "\n")
            n_kept += 1
    return n_total, n_removed, n_kept


def main():
    args = parse_args()

    if not args.output_vcf_gz.endswith(".vcf.gz"):
        sys.exit("[ERROR] --output-vcf-gz must end with .vcf.gz")

    # Step 1: build the lookup set from the filter VCF
    sys.stderr.write(f"[INFO] Loading filter variants from: {args.filter_vcf}\n")
    filter_keys, n_filter_records = load_filter_keys(args.filter_vcf)
    sys.stderr.write(f"[INFO] Filter VCF records loaded: {n_filter_records}\n")
    sys.stderr.write(f"[INFO] Filter allele keys loaded:  {len(filter_keys)}\n")

    # Step 2: filter the input into a temporary uncompressed VCF, then
    # bgzip + index it (the temp directory is removed automatically)
    with tempfile.TemporaryDirectory(prefix="vcf_overlap_filter_") as tmpdir:
        temp_vcf = os.path.join(tmpdir, "filtered.tmp.vcf")

        sys.stderr.write(f"[INFO] Filtering input VCF: {args.input_vcf}\n")
        n_total, n_removed, n_kept = filter_vcf_to_temp(
            input_vcf=args.input_vcf,
            filter_keys=filter_keys,
            temp_vcf=temp_vcf,
        )
        sys.stderr.write(f"[INFO] Total input records: {n_total}\n")
        sys.stderr.write(f"[INFO] Removed records:     {n_removed}\n")
        sys.stderr.write(f"[INFO] Kept records:        {n_kept}\n")

        sys.stderr.write(f"[INFO] Writing bgzipped output: {args.output_vcf_gz}\n")

        with open(args.output_vcf_gz, "wb") as fout:
            subprocess.check_call(["bgzip", "-c", temp_vcf], stdout=fout)
        subprocess.check_call(["tabix", "-f", "-p", "vcf", args.output_vcf_gz])

    sys.stderr.write(f"[DONE] Output VCF:   {args.output_vcf_gz}\n")
    sys.stderr.write(f"[DONE] Output index: {args.output_vcf_gz}.tbi\n")


if __name__ == "__main__":
    main()
