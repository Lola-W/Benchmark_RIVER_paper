#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Usage:
# python post_hoc_01.merge_benchmark_with_pileup_hg38_hg19.py \
#   --benchmark benchmark_snv_matrix.annot.with_hg19.tsv \
#   --hg38_root baseQ13_hg38/result \
#   --hg19_root baseQ13_hg19/result \
#   --out benchmark_snv_matrix.annot.with_hg19.plus_pileup.tsv \
#   --bulk_methods ILLUMINA,PPMSEQ,UDSEQ,NANOSEQ \
#   --ref_count_min 3 \
#   --upper_ci_thr 0.52

import argparse
import re
from pathlib import Path
from typing import List, Set

import numpy as np
import pandas as pd


METRIC_COLS = ["REF_COUNT", "ALT_COUNT", "MAF", "LOWER_CI", "UPPER_CI", "IS_GREATER"]


def norm_prefix(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", name.strip())
    return re.sub(r"_+", "_", s).strip("_").upper()


def norm_chr(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"^chr", "", x, flags=re.IGNORECASE)
    return x.upper()


def allowed_hg19_chrom(x: str) -> bool:
    x = norm_chr(x)
    if x == "X":
        return True
    if x.isdigit():
        return 1 <= int(x) <= 22
    return False


def make_key(df: pd.DataFrame, c, p, r, a) -> pd.Series:
    return (
        df[c].astype(str)
        .str.cat(df[p].astype(str), sep="|")
        .str.cat(df[r].astype(str), sep="|")
        .str.cat(df[a].astype(str), sep="|")
    )


def list_metric_files(root: Path, fname: str) -> List[Path]:
    files = []
    for d in sorted([p for p in root.iterdir() if p.is_dir()]):
        f = d / fname
        if f.exists() and f.stat().st_size > 0:
            files.append(f)
    return files


def read_metrics(path: Path, chrom_col="#CHROM", hg19=False) -> pd.DataFrame:
    usecols = [chrom_col, "POS", "REF", "ALT"] + METRIC_COLS
    df = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str)

    if hg19:
        df[chrom_col] = df[chrom_col].map(norm_chr)
        df = df[df[chrom_col].map(allowed_hg19_chrom)]

    df["__key__"] = make_key(df, chrom_col, "POS", "REF", "ALT")
    df = df[["__key__"] + METRIC_COLS]
    df = df.drop_duplicates("__key__", keep="first")
    return df


def merge_left(base: pd.DataFrame, dfm: pd.DataFrame, key: str, prefix: str) -> pd.DataFrame:
    dfm = dfm.rename(columns={c: f"{prefix}_{c}" for c in METRIC_COLS})
    out = base.merge(dfm, how="left", left_on=key, right_on="__key__")
    return out.drop(columns="__key__")


def normalize_missing(df: pd.DataFrame):
    for c in df.columns:
        if c.endswith("_REF_COUNT") or c.endswith("_ALT_COUNT"):
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
        elif c.endswith(("_MAF", "_LOWER_CI", "_UPPER_CI")):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        elif c.endswith("_IS_GREATER"):
            df[c] = df[c].where(df[c].notna(), np.nan)


def add_germline_flags(df, bulk_methods, ref_min, upper_ci_thr):
    flags = []

    for b in bulk_methods:
        refc = f"{b}_REF_COUNT"
        uci = f"{b}_UPPER_CI"
        flag = f"is_germline_{b}_UPPER_CI"

        if refc not in df.columns or uci not in df.columns:
            df[flag] = 0
        else:
            df[flag] = (
                (df[refc] >= ref_min)
                & (df[uci].notna())
                & (df[uci] > upper_ci_thr)
            ).astype(int)

        flags.append(flag)

    df["is_germline_any_bulk_upper_ci"] = (df[flags].sum(axis=1) > 0).astype(int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--hg38_root", required=True)
    ap.add_argument("--hg19_root", required=True)
    ap.add_argument("--hg38_filename", default="baseQ13.big_with_downsampling.maf_ci.tsv")
    ap.add_argument("--hg19_filename", default="baseQ13.big_with_downsampling.hg19.maf_ci.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bulk_methods", default="ILLUMINA,PPMSEQ,UDSEQ,NANOSEQ")
    ap.add_argument("--ref_count_min", type=int, default=3)
    ap.add_argument("--upper_ci_thr", type=float, default=0.52)
    args = ap.parse_args()

    bulk_methods = [b.strip().upper() for b in args.bulk_methods.split(",")]

    print("[INIT] loading benchmark", flush=True)
    bench = pd.read_csv(args.benchmark, sep="\t", dtype=str)
    print(f"[INIT] benchmark rows = {len(bench)}", flush=True)

    bench["__key_hg38__"] = make_key(bench, "chr", "pos", "ref", "alt")

    bench["hg19_chr_norm"] = bench["hg19_chr"].map(norm_chr)
    valid19 = (bench["hg19_unmapped"] == "0") & (bench["hg19_chr"] != ".")
    bench["__key_hg19__"] = ""
    bench.loc[valid19, "__key_hg19__"] = make_key(
        bench.loc[valid19],
        "hg19_chr_norm", "hg19_pos", "hg19_ref", "hg19_alt",
    )

    # -----------------------------
    # HG19 pileup presence flag
    # -----------------------------
    hg19_files_for_presence = list_metric_files(Path(args.hg19_root), args.hg19_filename)
    print(f"[HG19_PRESENCE] total files = {len(hg19_files_for_presence)}", flush=True)

    hg19_presence: Set[str] = set()
    for i, f in enumerate(hg19_files_for_presence, 1):
        src = f.parent.name
        print(f"[HG19_PRESENCE] ({i}/{len(hg19_files_for_presence)}) reading {src}", flush=True)

        dfm = read_metrics(f, chrom_col="#CHROM", hg19=True)
        hg19_presence.update(dfm["__key__"].tolist())

    print(f"[HG19_PRESENCE] unique keys = {len(hg19_presence)}", flush=True)

    bench["hg19_pileup_unmapped"] = "1"
    present = valid19 & bench["__key_hg19__"].isin(hg19_presence)
    bench.loc[present, "hg19_pileup_unmapped"] = "0"


    # HG38
    hg38_files = list_metric_files(Path(args.hg38_root), args.hg38_filename)
    print(f"[HG38] total files = {len(hg38_files)}", flush=True)

    for i, f in enumerate(hg38_files, 1):
        src = f.parent.name
        print(f"[HG38] ({i}/{len(hg38_files)}) merging {src}", flush=True)
        prefix = norm_prefix(src)
        dfm = read_metrics(f, "#CHROM", hg19=False)
        bench = merge_left(bench, dfm, "__key_hg38__", prefix)

    # HG19
    hg19_files = list_metric_files(Path(args.hg19_root), args.hg19_filename)
    print(f"[HG19] total files = {len(hg19_files)}", flush=True)

    for i, f in enumerate(hg19_files, 1):
        src = f.parent.name
        print(f"[HG19] ({i}/{len(hg19_files)}) merging {src}", flush=True)
        prefix = "HG19_" + norm_prefix(src)
        dfm = read_metrics(f, "#CHROM", hg19=True)
        bench = merge_left(bench, dfm, "__key_hg19__", prefix)

    print("[POST] normalizing missing values", flush=True)
    bench.drop(columns=["__key_hg38__", "__key_hg19__", "hg19_chr_norm"], inplace=True)
    normalize_missing(bench)

    print("[POST] computing germline flags", flush=True)
    add_germline_flags(
        bench,
        bulk_methods,
        args.ref_count_min,
        args.upper_ci_thr,
    )

    print("[WRITE] writing output", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    bench.to_csv(args.out, sep="\t", index=False)

    print("[DONE]", flush=True)


if __name__ == "__main__":
    main()

