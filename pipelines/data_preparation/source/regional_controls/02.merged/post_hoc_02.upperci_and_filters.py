#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Usage:
# python post_hoc_02.upperci_and_filters.py \
#   --in benchmark_snv_matrix.annot.with_hg19.plus_pileup.tsv \
#   --out_prefix benchmark_snv_matrix.annot.with_hg19 \
#   --vcf_outdir sample_naive_filtered_vcf \
#   --write_primary_table \
#   --make_anybulk_table \
#   --apply_primary_filters \
#   --bulk_methods ILLUMINA,PPMSEQ,UDSEQ,NANOSEQ
#
# Notes:
# - Primary filters = hg19_pileup_unmapped==0 and (blacklist_any, gnomad_common, in_known_indel_any, is_germline)==0
# - If --write_primary_table is set, we write:
#     <out_prefix>.primary_filter.tsv
# - If --make_anybulk_table is set, we write:
#     <out_prefix>.any_bulk_upper_ci_filter.tsv  (applies is_germline_any_bulk_upper_ci==0 on top of primary-filtered rows)
# - Per-caller VCFs are always written from PRIMARY-filtered rows:
#     bulk callers additionally apply caller-specific sample-naive upperCI flags

import argparse
from pathlib import Path
import pandas as pd

SAMPLE_NAIVE_FLAG_BY_CALLER = {
    "Illumina": "is_germline_ILLUMINA_UPPER_CI",
    "ppmSeq_multiread": "is_germline_PPMSEQ_UPPER_CI",
    "ppmSeq_putative_multiread": "is_germline_PPMSEQ_UPPER_CI",
    "ppmSeq_singleton_HC": "is_germline_PPMSEQ_UPPER_CI",
    "UDSeq": "is_germline_UDSEQ_UPPER_CI",
    "NanoSeq": "is_germline_NANOSEQ_UPPER_CI",
}

CALLERS_TO_EXPORT = [
    "Illumina",
    "ppmSeq_multiread",
    "ppmSeq_putative_multiread",
    "ppmSeq_singleton_HC",
    "UDSeq",
    "NanoSeq",
    "HiDEFseq_double_stranded",
    "HiDEFseq_single_stranded",
    "PTA_single_neuron_01",
    "PTA_single_neuron_02",
    "PTA_single_neuron_03",
    "PTA_single_neuron_04",
    "PTA_single_neuron_05",
    "PTA_single_neuron_06",
    "PTA_single_neuron_07",
    "PTA_single_neuron_08",
    "PTA_single_neuron_09",
    "PTA_single_neuron_10",
    "PTA_single_neuron_11",
    "PTA_single_neuron_12",
    "PTA_single_neuron_13",
    "PTA_single_neuron_14",
    "PTA_single_neuron_15",
    "PTA_single_neuron_16",
    "AmpliSeq",
]

PRIMARY_FILTER_COLS = [
    "hg19_pileup_unmapped",
    "blacklist_any",
    "gnomad_common",
    "in_known_indel_any",
    "is_germline",
]


def ensure_cols(df, cols):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")


def apply_primary_filters(df: pd.DataFrame) -> pd.DataFrame:
    # Enforce:
    # hg19_pileup_unmapped == 0
    # blacklist_any == 0
    # gnomad_common == 0
    # in_known_indel_any == 0
    # is_germline == 0
    m = (
        (df["hg19_pileup_unmapped"].astype(str) == "0") &
        (df["blacklist_any"].astype(str) == "0") &
        (df["gnomad_common"].astype(str) == "0") &
        (df["in_known_indel_any"].astype(str) == "0") &
        (df["is_germline"].astype(str) == "0")
    )
    return df[m].copy()


def write_vcf(df: pd.DataFrame, out_vcf: Path, source_tag: str):
    out_vcf.parent.mkdir(parents=True, exist_ok=True)
    with out_vcf.open("w") as w:
        w.write("##fileformat=VCFv4.2\n")
        w.write(f"##source={source_tag}\n")
        w.write('##INFO=<ID=SRC,Number=1,Type=String,Description="Caller name">\n')
        w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, r in df.iterrows():
            chrom = str(r["chr"])
            pos = str(r["pos"])
            ref = str(r["ref"])
            alt = str(r["alt"])
            info = f"SRC={source_tag}"
            w.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV (plus_pileup)")
    ap.add_argument("--out_prefix", required=True, help="Prefix for output files")
    ap.add_argument("--vcf_outdir", default="sample_naive_filtered_vcf", help="Output dir for per-caller VCFs")
    ap.add_argument("--apply_primary_filters", action="store_true",
                    help="Apply primary filters before writing outputs (recommended)")
    ap.add_argument("--write_primary_table", action="store_true",
                    help="Write <out_prefix>.primary_filter.tsv (primary-filtered table)")
    ap.add_argument("--make_anybulk_table", action="store_true",
                    help="Write <out_prefix>.any_bulk_upper_ci_filter.tsv (primary + is_germline_any_bulk_upper_ci==0)")
    ap.add_argument("--bulk_methods", default="ILLUMINA,PPMSEQ,UDSEQ,NANOSEQ",
                    help="Kept for compatibility; not required for current logic")
    args = ap.parse_args()

    print(f"[LOAD] {args.inp}", flush=True)
    df = pd.read_csv(args.inp, sep="\t", dtype=str)
    print(f"[LOAD] rows = {len(df)}", flush=True)

    ensure_cols(df, ["chr", "pos", "ref", "alt"])
    ensure_cols(df, PRIMARY_FILTER_COLS)

    # Primary filter
    if args.apply_primary_filters:
        n_before = len(df)
        df_primary = apply_primary_filters(df)
        print(f"[FILTER] primary filters applied: {n_before} -> {len(df_primary)}", flush=True)
    else:
        # Still allow writing without filtering if user insists
        df_primary = df.copy()
        print("[FILTER] primary filters NOT applied (as requested)", flush=True)

    # Write primary-filtered table
    if args.write_primary_table:
        out_primary = f"{args.out_prefix}.primary_filter.tsv"
        df_primary.to_csv(out_primary, sep="\t", index=False)
        print(f"[WRITE] {out_primary} rows={len(df_primary)}", flush=True)

    # Write any-bulk filtered table (on top of primary)
    if args.make_anybulk_table:
        if "is_germline_any_bulk_upper_ci" not in df_primary.columns:
            raise RuntimeError("Missing column: is_germline_any_bulk_upper_ci (needed for any-bulk table)")
        n_before = len(df_primary)
        df_anybulk = df_primary[df_primary["is_germline_any_bulk_upper_ci"].astype(str) != "1"].copy()
        out_any = f"{args.out_prefix}.any_bulk_upper_ci_filter.tsv"
        df_anybulk.to_csv(out_any, sep="\t", index=False)
        print(f"[WRITE] {out_any} rows={len(df_anybulk)} (from {n_before})", flush=True)

    # Export per-caller VCFs from PRIMARY-filtered rows
    outdir = Path(args.vcf_outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[VCF] exporting VCFs to {outdir}", flush=True)

    for i, caller in enumerate(CALLERS_TO_EXPORT, 1):
        if caller not in df_primary.columns:
            print(f"[VCF] ({i}/{len(CALLERS_TO_EXPORT)}) skip {caller} (missing column)", flush=True)
            continue

        sub = df_primary[df_primary[caller].astype(str) == "1"].copy()

        flag = SAMPLE_NAIVE_FLAG_BY_CALLER.get(caller, None)
        if flag is not None:
            if flag in sub.columns:
                n0 = len(sub)
                sub = sub[sub[flag].astype(str) != "1"].copy()
                print(f"[VCF] ({i}/{len(CALLERS_TO_EXPORT)}) {caller}: call=1 {n0} -> upperCI({flag}) -> {len(sub)}", flush=True)
            else:
                print(f"[VCF] ({i}/{len(CALLERS_TO_EXPORT)}) {caller}: flag {flag} missing; no upperCI filter applied", flush=True)
        else:
            print(f"[VCF] ({i}/{len(CALLERS_TO_EXPORT)}) {caller}: no sample-naive upperCI; primary-only rows={len(sub)}", flush=True)

        out_vcf = outdir / f"{caller}.sample_naive_filtered.vcf"
        write_vcf(sub, out_vcf, source_tag=caller)

    print("[DONE]", flush=True)


if __name__ == "__main__":
    main()

