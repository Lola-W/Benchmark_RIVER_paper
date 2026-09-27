#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# reorder_add_flags.py
#
# Usage:
# python reorder_add_flags.py \
#   --in_tsv  benchmark_with_downsampling_quant_added.reordered.tsv \
#   --out_tsv benchmark_with_downsampling_quant_added.flags.tsv \
#   --alt_min 1 \
#   --depth_min 3 \
#   --uci_thr 0.52 \
#   --uci_lab 052 \
#   --drop_old_germline_flags

import argparse
import re
from pathlib import Path
from typing import List, Set, Tuple

import numpy as np
import pandas as pd

METRIC_COLS = ["REF_COUNT", "ALT_COUNT", "MAF", "LOWER_CI", "UPPER_CI", "IS_GREATER"]

def find_metric_prefixes(cols: List[str]) -> Set[str]:
    prefixes = set()
    pat = re.compile(r"^(.*)_(REF_COUNT|ALT_COUNT|MAF|LOWER_CI|UPPER_CI|IS_GREATER)$")
    for c in cols:
        m = pat.match(c)
        if m:
            prefixes.add(m.group(1))
    return prefixes

def ensure_missing_rule(df: pd.DataFrame, prefix: str) -> None:
    for m in METRIC_COLS:
        col = f"{prefix}_{m}"
        if col not in df.columns:
            df[col] = 0 if m in ("REF_COUNT", "ALT_COUNT") else pd.NA

    refc = f"{prefix}_REF_COUNT"
    altc = f"{prefix}_ALT_COUNT"
    df[refc] = pd.to_numeric(df[refc], errors="coerce").fillna(0).astype(int)
    df[altc] = pd.to_numeric(df[altc], errors="coerce").fillna(0).astype(int)

    for m in ("MAF", "LOWER_CI", "UPPER_CI"):
        col = f"{prefix}_{m}"
        df[col] = pd.to_numeric(df[col], errors="coerce")

    isg = f"{prefix}_IS_GREATER"
    df[isg] = df[isg].replace({"": pd.NA})

def method_order_key(prefix: str) -> Tuple[int, int, str]:
    base_rank = {"ILLUMINA": 1, "PPMSEQ": 2, "UDSEQ": 3, "NANOSEQ": 4, "HIDEFSEQ": 5}
    m = re.match(r"^(.*)_(\d+)x$", prefix)
    if m:
        base = m.group(1)
        depth = int(m.group(2))
        depth_rank = {30: 1, 100: 2, 200: 3}.get(depth, 9)
        br = base_rank.get(base, 99)
        return (br, depth_rank, prefix)
    br = base_rank.get(prefix, 99)
    return (br, 0, prefix)

def make_flag_cols_for_prefix(df: pd.DataFrame, prefix: str, alt_min: int, depth_min: int, uci_thr: float, uci_lab: str) -> Tuple[str, str]:
    ensure_missing_rule(df, prefix)
    refc = df[f"{prefix}_REF_COUNT"]
    altc = df[f"{prefix}_ALT_COUNT"]
    uci  = df[f"{prefix}_UPPER_CI"]

    col_uci = f"is_maf_upper_ci_ge_{uci_lab}_{prefix}"
    cond = (altc >= alt_min) & ((altc + refc) >= depth_min) & (uci.notna()) & (uci > uci_thr)
    df[col_uci] = cond.astype(int)

    col_alt = f"is_alt_count_ge_{alt_min}_{prefix}"
    df[col_alt] = (altc >= alt_min).astype(int)
    return col_uci, col_alt

def pta_summary_maf(df: pd.DataFrame, alt_min: int, depth_min: int, uci_thr: float, uci_lab: str, pta_prefixes: List[str]) -> Tuple[str, str, str]:
    for p in pta_prefixes:
        ensure_missing_rule(df, p)

    count = np.zeros(len(df), dtype=int)
    for p in pta_prefixes:
        refc = df[f"{p}_REF_COUNT"].to_numpy()
        altc = df[f"{p}_ALT_COUNT"].to_numpy()
        uci  = df[f"{p}_UPPER_CI"].to_numpy()
        ok = (altc >= alt_min) & ((altc + refc) >= depth_min) & np.isfinite(uci) & (uci > uci_thr)
        count += ok.astype(int)

    col_num  = f"is_maf_upper_ci_ge_{uci_lab}_PTA_number"
    col_sing = f"is_maf_upper_ci_ge_{uci_lab}_PTA_singleton"
    col_mult = f"is_maf_upper_ci_ge_{uci_lab}_PTA_multiple"
    df[col_num]  = count.astype(int)
    df[col_sing] = (df[col_num] == 1).astype(int)
    df[col_mult] = (df[col_num] >= 2).astype(int)
    return col_num, col_sing, col_mult

def pta_summary_alt(df: pd.DataFrame, alt_min: int, pta_prefixes: List[str]) -> Tuple[str, str, str]:
    for p in pta_prefixes:
        ensure_missing_rule(df, p)
    count = np.zeros(len(df), dtype=int)
    for p in pta_prefixes:
        altc = df[f"{p}_ALT_COUNT"].to_numpy()
        count += (altc >= alt_min).astype(int)

    col_num  = f"is_alt_count_ge_{alt_min}_PTA_number"
    col_sing = f"is_alt_count_ge_{alt_min}_PTA_singleton"
    col_mult = f"is_alt_count_ge_{alt_min}_PTA_multiple"
    df[col_num]  = count.astype(int)
    df[col_sing] = (df[col_num] == 1).astype(int)
    df[col_mult] = (df[col_num] >= 2).astype(int)
    return col_num, col_sing, col_mult

def other_brain_region_flags(df: pd.DataFrame, alt_min: int) -> Tuple[str, str]:
    targets = [
        "HG19_ILLUMINA_B_L_F_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_O_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_PF_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_P_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_F_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_O_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_PF_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_P_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_T_LRG_ALT_COUNT",
    ]
    present = [c for c in targets if c in df.columns]
    col_flag  = f"is_alt_count_ge_{alt_min}_other_brain_region"
    col_which = f"is_alt_count_ge_{alt_min}_other_brain_region_which"
    if not present:
        df[col_flag] = 0
        df[col_which] = ""
        return col_flag, col_which

    for c in present:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    hits = (df[present] >= alt_min)
    which = []
    for _, row in hits.iterrows():
        cols = [c for c in present if row[c]]
        tags = [c.replace("HG19_ILLUMINA_", "").replace("_LRG_ALT_COUNT", "") for c in cols]
        which.append(",".join(tags))

    df[col_flag] = hits.any(axis=1).astype(int)
    df[col_which] = which
    return col_flag, col_which

def ctrl_platform_from_tag(tag: str) -> str:
    if tag.startswith("HG19_ILLUMINA_"):
        return "ILLUMINA"
    if tag.startswith("HG19_PTA_"):
        return "PTA"
    for p in ["ILLUMINA", "PPMSEQ", "UDSEQ", "NANOSEQ", "HIDEFSEQ"]:
        if tag.startswith(p + "_"):
            return p
    return tag.split("_", 1)[0]

def unrelated_control_flags(df: pd.DataFrame, alt_min: int) -> Tuple[str, str, str]:
    alt_cols = [c for c in df.columns if c.endswith("_ALT_COUNT")]

    def is_ctrl(c: str) -> bool:
        if c.startswith("HG19_PTA_") and c.endswith("_ALT_COUNT"):
            return True
        if re.search(r"(6566|7669|8420)", c):
            return True
        return False

    ctrl_cols = [c for c in alt_cols if is_ctrl(c)]

    col_flag = f"is_alt_count_ge_{alt_min}_unrelated_control"
    col_plat = f"is_alt_count_ge_{alt_min}_unrelated_control_platform"
    col_which = f"is_alt_count_ge_{alt_min}_unrelated_control_which"

    if not ctrl_cols:
        df[col_flag] = 0
        df[col_plat] = ""
        df[col_which] = ""
        return col_flag, col_plat, col_which

    for c in ctrl_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    hits = (df[ctrl_cols] >= alt_min)

    which_list = []
    plat_list = []
    for _, row in hits.iterrows():
        cols = [c for c in ctrl_cols if row[c]]
        tags = [c.replace("_ALT_COUNT", "") for c in cols]
        which_list.append(",".join(tags))
        plats = sorted({ctrl_platform_from_tag(t) for t in tags})
        plat_list.append(",".join(plats))

    df[col_flag] = hits.any(axis=1).astype(int)
    df[col_plat] = plat_list
    df[col_which] = which_list

    df.loc[df[col_flag] == 0, col_plat] = ""
    df.loc[df[col_flag] == 0, col_which] = ""
    return col_flag, col_plat, col_which

def reorder_new_columns(df: pd.DataFrame, uci_lab: str, alt_min: int, bench_prefixes_sorted: List[str], has_pta: bool) -> List[str]:
    maf_cols = [f"is_maf_upper_ci_ge_{uci_lab}_{p}" for p in bench_prefixes_sorted]
    alt_cols = [f"is_alt_count_ge_{alt_min}_{p}" for p in bench_prefixes_sorted]
    ordered = []
    ordered.extend([c for c in maf_cols if c in df.columns])
    if has_pta:
        ordered.extend([c for c in [
            f"is_maf_upper_ci_ge_{uci_lab}_PTA_number",
            f"is_maf_upper_ci_ge_{uci_lab}_PTA_singleton",
            f"is_maf_upper_ci_ge_{uci_lab}_PTA_multiple",
        ] if c in df.columns])
    ordered.extend([c for c in alt_cols if c in df.columns])
    if has_pta:
        ordered.extend([c for c in [
            f"is_alt_count_ge_{alt_min}_PTA_number",
            f"is_alt_count_ge_{alt_min}_PTA_singleton",
            f"is_alt_count_ge_{alt_min}_PTA_multiple",
        ] if c in df.columns])
    ordered.extend([c for c in [
        f"is_alt_count_ge_{alt_min}_other_brain_region",
        f"is_alt_count_ge_{alt_min}_other_brain_region_which",
        f"is_alt_count_ge_{alt_min}_unrelated_control",
        f"is_alt_count_ge_{alt_min}_unrelated_control_platform",
        f"is_alt_count_ge_{alt_min}_unrelated_control_which",
    ] if c in df.columns])
    return ordered

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--alt_min", type=int, default=1)
    ap.add_argument("--depth_min", type=int, default=3)
    ap.add_argument("--uci_thr", type=float, default=0.52)
    ap.add_argument("--uci_lab", required=True)  # <-- FORCE LABEL FROM CLI
    ap.add_argument("--drop_old_germline_flags", action="store_true")
    args = ap.parse_args()

    print(f"[LOAD] {args.in_tsv}", flush=True)
    df = pd.read_csv(args.in_tsv, sep="\t", dtype=str, low_memory=False)
    print(f"[LOAD] rows={len(df):,} cols={len(df.columns):,}", flush=True)
    print(f"[PARAM] alt_min={args.alt_min} depth_min={args.depth_min} uci_thr={args.uci_thr} uci_lab={args.uci_lab}", flush=True)

    if args.drop_old_germline_flags:
        drop_cols = [
            "is_germline_ILLUMINA_UPPER_CI",
            "is_germline_PPMSEQ_UPPER_CI",
            "is_germline_UDSEQ_UPPER_CI",
            "is_germline_NANOSEQ_UPPER_CI",
            # also handle accidental duplicated columns from previous merges
            "is_germline_ILLUMINA_UPPER_CI.1",
            "is_germline_PPMSEQ_UPPER_CI.1",
            "is_germline_UDSEQ_UPPER_CI.1",
            "is_germline_NANOSEQ_UPPER_CI.1",
        ]
        exist = [c for c in drop_cols if c in df.columns]
        if exist:
            df.drop(columns=exist, inplace=True)
            print(f"[DROP] removed {len(exist)} old germline flag cols", flush=True)

    prefixes = sorted(find_metric_prefixes(list(df.columns)))

    bench_prefixes = []
    for p in prefixes:
        if p.startswith("HG19_"):
            continue
        if p.startswith("PTA_SINGLE_NEURON_"):
            continue
        if re.search(r"(6566|7669|8420)", p):
            continue
        bench_prefixes.append(p)

    bench_prefixes_sorted = sorted(bench_prefixes, key=method_order_key)
    print(f"[BENCH] metric prefixes detected={len(prefixes)} ; benchmark targets={len(bench_prefixes_sorted)}", flush=True)

    # create bench flags
    for p in bench_prefixes_sorted:
        make_flag_cols_for_prefix(df, p, args.alt_min, args.depth_min, args.uci_thr, args.uci_lab)

    # PTA prefixes present?
    pta_prefixes = [f"PTA_SINGLE_NEURON_{i:02d}" for i in range(1, 17)]
    pta_prefixes = [p for p in pta_prefixes if (f"{p}_ALT_COUNT" in df.columns or f"{p}_REF_COUNT" in df.columns or f"{p}_UPPER_CI" in df.columns)]
    has_pta = len(pta_prefixes) > 0

    if has_pta:
        pta_summary_maf(df, args.alt_min, args.depth_min, args.uci_thr, args.uci_lab, pta_prefixes)
        pta_summary_alt(df, args.alt_min, pta_prefixes)
        print(f"[PTA] prefixes used={len(pta_prefixes)}", flush=True)
    else:
        print("[PTA] no PTA prefixes found, skipped", flush=True)

    other_brain_region_flags(df, args.alt_min)
    unrelated_control_flags(df, args.alt_min)

    # reorder columns: original first + ordered new flags
    new_flag_cols = [c for c in df.columns if c.startswith("is_maf_upper_ci_ge_") or c.startswith("is_alt_count_ge_")]
    ordered_new = reorder_new_columns(df, args.uci_lab, args.alt_min, bench_prefixes_sorted, has_pta)

    original_cols = [c for c in df.columns if c not in set(new_flag_cols)]
    extra_new = [c for c in new_flag_cols if c not in set(ordered_new)]
    df = df[original_cols + ordered_new + extra_new]

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    print(f"[WRITE] {args.out_tsv}", flush=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[DONE] rows={len(df):,} cols={len(df.columns):,} ordered_new={len(ordered_new)} extra_new={len(extra_new)}", flush=True)

if __name__ == "__main__":
    main()

