#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path
from typing import List, Tuple, Set

import numpy as np
import pandas as pd


METRIC_COLS = ["REF_COUNT", "ALT_COUNT", "MAF", "LOWER_CI", "UPPER_CI", "IS_GREATER"]


def thr_label(x: float) -> str:
    """
    Convert threshold to 3-digit label (percent*100, zero-padded).
      0.52 -> "052"
      0.64 -> "064"
      0.7  -> "070"
      1.0  -> "100"
    """
    lab = int(round(x * 100))
    if lab < 0:
        lab = 0
    if lab > 999:
        lab = 999
    return f"{lab:03d}"


def find_metric_prefixes(cols: List[str]) -> Set[str]:
    """
    Return prefixes that have at least one of the metric cols.
    Prefix definition: <PREFIX>_<METRIC>, where METRIC in METRIC_COLS.
    """
    prefixes = set()
    pat = re.compile(r"^(.*)_(REF_COUNT|ALT_COUNT|MAF|LOWER_CI|UPPER_CI|IS_GREATER)$")
    for c in cols:
        m = pat.match(c)
        if m:
            prefixes.add(m.group(1))
    return prefixes


def ensure_missing_rule(df: pd.DataFrame, prefix: str) -> None:
    """
    Missing rule for a prefix:
      - REF_COUNT, ALT_COUNT: missing -> 0
      - MAF, LOWER_CI, UPPER_CI, IS_GREATER: missing -> NA

    Also coerces:
      - *_REF_COUNT, *_ALT_COUNT -> int (NA -> 0)
      - *_MAF, *_LOWER_CI, *_UPPER_CI -> float (NaN allowed)
      - *_IS_GREATER -> object (NA allowed)
    """
    # Create missing columns
    for m in METRIC_COLS:
        col = f"{prefix}_{m}"
        if col not in df.columns:
            if m in ("REF_COUNT", "ALT_COUNT"):
                df[col] = 0
            else:
                df[col] = pd.NA

    # Coerce types
    refc = f"{prefix}_REF_COUNT"
    altc = f"{prefix}_ALT_COUNT"
    df[refc] = pd.to_numeric(df[refc], errors="coerce").fillna(0).astype(int)
    df[altc] = pd.to_numeric(df[altc], errors="coerce").fillna(0).astype(int)

    for m in ("MAF", "LOWER_CI", "UPPER_CI"):
        col = f"{prefix}_{m}"
        df[col] = pd.to_numeric(df[col], errors="coerce")  # keep NaN

    isg = f"{prefix}_IS_GREATER"
    df[isg] = df[isg].replace({"": pd.NA})


def make_flag_cols_for_prefix(
    df: pd.DataFrame,
    prefix: str,
    alt_min: int,
    depth_min: int,
    uci_thr: float,
    uci_lab: str,
) -> Tuple[str, str]:
    """
    Adds:
      is_alt_count_ge_<alt_min>_<PREFIX>
      is_maf_upper_ci_ge_<uci_lab>_<PREFIX>
    """
    ensure_missing_rule(df, prefix)

    refc = df[f"{prefix}_REF_COUNT"]
    altc = df[f"{prefix}_ALT_COUNT"]
    uci = df[f"{prefix}_UPPER_CI"]

    col_alt = f"is_alt_count_ge_{alt_min}_{prefix}"
    df[col_alt] = (altc >= alt_min).astype(int)

    # Condition: ALT>=alt_min AND ALT+REF>=depth_min AND UPPER_CI>uci_thr
    cond = (altc >= alt_min) & ((altc + refc) >= depth_min) & (uci.notna()) & (uci > uci_thr)
    col_uci = f"is_maf_upper_ci_ge_{uci_lab}_{prefix}"
    df[col_uci] = cond.astype(int)

    return col_alt, col_uci


def pta_summary(
    df: pd.DataFrame,
    alt_min: int,
    depth_min: int,
    uci_thr: float,
    uci_lab: str,
    pta_prefixes: List[str],
) -> Tuple[str, str, str]:
    """
    PTA prefixes에 대해 조건 만족하는 개수를 세어 number/singleton/multiple 생성.
    """
    for p in pta_prefixes:
        ensure_missing_rule(df, p)

    count = np.zeros(len(df), dtype=int)

    for p in pta_prefixes:
        refc = df[f"{p}_REF_COUNT"].to_numpy()
        altc = df[f"{p}_ALT_COUNT"].to_numpy()
        uci = df[f"{p}_UPPER_CI"].to_numpy()

        # uci is float array with NaN; np.isfinite handles NaN
        ok = (altc >= alt_min) & ((altc + refc) >= depth_min) & np.isfinite(uci) & (uci > uci_thr)
        count += ok.astype(int)

    col_num = f"is_maf_upper_ci_ge_{uci_lab}_PTA_UPPER_CI_number"
    col_sing = f"is_maf_upper_ci_ge_{uci_lab}_PTA_UPPER_CI_singleton"
    col_mult = f"is_maf_upper_ci_ge_{uci_lab}_PTA_UPPER_CI_multiple"

    df[col_num] = count.astype(int)
    df[col_sing] = (df[col_num] == 1).astype(int)
    df[col_mult] = (df[col_num] >= 2).astype(int)
    return col_num, col_sing, col_mult


def other_brain_region_flags(df: pd.DataFrame, alt_min: int) -> Tuple[str, str]:
    """
    same sample different region:
      - 제외: B_L_T
      - 나머지에서 ALT>=alt_min 하나라도 있으면 flag=1, which=comma list (B_L_F,B_R_T,...)
    """
    targets = [
        "HG19_ILLUMINA_B_L_F_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_O_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_PF_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_L_P_LRG_ALT_COUNT",
        # "HG19_ILLUMINA_B_L_T_LRG_ALT_COUNT",  # excluded
        "HG19_ILLUMINA_B_R_F_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_O_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_PF_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_P_LRG_ALT_COUNT",
        "HG19_ILLUMINA_B_R_T_LRG_ALT_COUNT",
    ]

    present = [c for c in targets if c in df.columns]
    col_flag = f"is_alt_count_ge_{alt_min}_other_brain_region"
    col_which = f"is_alt_count_ge_{alt_min}_other_brain_region_which"

    if not present:
        df[col_flag] = 0
        df[col_which] = ""
        return col_flag, col_which

    for c in present:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    hits = (df[present] >= alt_min)

    # vectorized-ish string build
    which = []
    for _, row in hits.iterrows():
        cols = [c for c in present if row[c]]
        tags = [c.replace("HG19_ILLUMINA_", "").replace("_LRG_ALT_COUNT", "") for c in cols]
        which.append(",".join(tags))

    df[col_flag] = hits.any(axis=1).astype(int)
    df[col_which] = which
    return col_flag, col_which


def unrelated_control_flags(df: pd.DataFrame, alt_min: int) -> Tuple[str, str]:
    """
    different sample negative control:
      - ALT_COUNT columns whose name contains 6566 or 7669 or 8420
      - OR starts with HG19_PTA_
    """
    alt_cols = [c for c in df.columns if c.endswith("_ALT_COUNT")]

    def is_ctrl(c: str) -> bool:
        if c.startswith("HG19_PTA_"):
            return True
        if re.search(r"(6566|7669|8420)", c):
            return True
        return False

    ctrl_cols = [c for c in alt_cols if is_ctrl(c)]

    col_flag = f"is_alt_count_ge_{alt_min}_unrelated_control"
    col_which = f"is_alt_count_ge_{alt_min}_unrelated_control_which"

    if not ctrl_cols:
        df[col_flag] = 0
        df[col_which] = ""
        return col_flag, col_which

    for c in ctrl_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)

    hits = (df[ctrl_cols] >= alt_min)

    which = []
    for _, row in hits.iterrows():
        cols = [c for c in ctrl_cols if row[c]]
        tags = [c.replace("_ALT_COUNT", "") for c in cols]
        which.append(",".join(tags))

    df[col_flag] = hits.any(axis=1).astype(int)
    df[col_which] = which
    return col_flag, col_which


def drop_old_germline_flags(df: pd.DataFrame) -> int:
    """
    Drop old germline flags, including pandas-mangled duplicates like '.1', '.2', ...
    """
    pat = re.compile(r"^is_germline_(ILLUMINA|PPMSEQ|UDSEQ|NANOSEQ)_UPPER_CI(\.\d+)?$")
    pat_any = re.compile(r"^is_germline_any_bulk_upper_ci(\.\d+)?$")

    to_drop = [c for c in df.columns if pat.match(c) or pat_any.match(c)]
    if to_drop:
        df.drop(columns=to_drop, inplace=True)
    return len(to_drop)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--alt_min", type=int, default=1)
    ap.add_argument("--depth_min", type=int, default=3)
    ap.add_argument("--uci_thr", type=float, default=0.52)
    ap.add_argument("--drop_old_germline_flags", action="store_true")
    args = ap.parse_args()

    uci_lab = thr_label(args.uci_thr)

    print(f"[LOAD] {args.in_tsv}", flush=True)
    df = pd.read_csv(args.in_tsv, sep="\t", dtype=str, low_memory=False)
    print(f"[LOAD] rows={len(df):,} cols={len(df.columns):,}", flush=True)

    if args.drop_old_germline_flags:
        n = drop_old_germline_flags(df)
        print(f"[DROP] removed {n} old germline flag cols (incl. .1/.2)", flush=True)

    # -----------------------------
    # SAME SAMPLE BENCHMARK prefixes (auto-detect)
    # exclude: HG19_* , PTA_SINGLE_NEURON_* (handled separately), and 6566/7669/8420 (controls)
    # -----------------------------
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

    print(f"[BENCH] metric prefixes detected={len(prefixes)} ; benchmark targets={len(bench_prefixes)}", flush=True)

    new_cols: List[str] = []

    for p in bench_prefixes:
        c1, c2 = make_flag_cols_for_prefix(
            df,
            p,
            alt_min=args.alt_min,
            depth_min=args.depth_min,
            uci_thr=args.uci_thr,
            uci_lab=uci_lab,
        )
        new_cols.extend([c1, c2])

    # -----------------------------
    # PTA summary (01-16)
    # -----------------------------
    pta_prefixes = [f"PTA_SINGLE_NEURON_{i:02d}" for i in range(1, 17)]
    pta_prefixes = [
        p for p in pta_prefixes
        if (f"{p}_ALT_COUNT" in df.columns) or (f"{p}_REF_COUNT" in df.columns) or (f"{p}_UPPER_CI" in df.columns)
    ]

    if pta_prefixes:
        cnum, cs, cm = pta_summary(
            df,
            alt_min=args.alt_min,
            depth_min=args.depth_min,
            uci_thr=args.uci_thr,
            uci_lab=uci_lab,
            pta_prefixes=pta_prefixes,
        )
        new_cols.extend([cnum, cs, cm])
        print(f"[PTA] prefixes used={len(pta_prefixes)}", flush=True)
    else:
        print("[PTA] no PTA prefixes found, skipped", flush=True)

    # -----------------------------
    # other brain region + unrelated controls
    # -----------------------------
    cflag, cwhich = other_brain_region_flags(df, alt_min=args.alt_min)
    new_cols.extend([cflag, cwhich])

    uflag, uwhich = unrelated_control_flags(df, alt_min=args.alt_min)
    new_cols.extend([uflag, uwhich])

    # -----------------------------
    # WRITE
    # -----------------------------
    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"[WRITE] {args.out_tsv}", flush=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[DONE] uci_lab={uci_lab} added_cols={len(new_cols)}", flush=True)


if __name__ == "__main__":
    main()

