#python fix_cols_and_fill.py \
#  --in_tsv postqc.slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.drop_allzero_call.tsv \
#  --out_tsv postqc.slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.drop_allzero_call.fixed.tsv \
#  --verbose

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path

import pandas as pd


METRIC_SUFFIX_ORDER = [
    "REF_COUNT",
    "ALT_COUNT",
    "MAF",
    "LOWER_CI",
    "UPPER_CI",
    "IS_GREATER",
]

METRIC_RE = re.compile(r"^(?P<prefix>.+)_(?P<suffix>REF_COUNT|ALT_COUNT|MAF|LOWER_CI|UPPER_CI|IS_GREATER)$")


def parse_args():
    p = argparse.ArgumentParser(
        description="Reorder pileup metric columns per prefix and fill empty fields with desired placeholders."
    )
    p.add_argument("--in_tsv", required=True, help="Input TSV")
    p.add_argument("--out_tsv", required=True, help="Output TSV")
    p.add_argument(
        "--verbose", action="store_true", help="Print summary to stderr"
    )
    return p.parse_args()


def is_empty(x) -> bool:
    # Treat NaN/None/empty/whitespace as empty
    if x is None:
        return True
    if isinstance(x, float) and pd.isna(x):
        return True
    s = str(x)
    return s.strip() == ""


def main():
    args = parse_args()

    in_path = Path(args.in_tsv)
    out_path = Path(args.out_tsv)

    # Read as strings to preserve exact tokens, then handle empty explicitly
    df = pd.read_csv(in_path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)

    cols = list(df.columns)

    # Identify metric columns and prefix order as they appear in the header
    prefix_order = []
    metric_cols_by_prefix = {}
    metric_indices = []

    for idx, c in enumerate(cols):
        m = METRIC_RE.match(c)
        if not m:
            continue
        metric_indices.append(idx)
        pref = m.group("prefix")
        suf = m.group("suffix")
        if pref not in metric_cols_by_prefix:
            metric_cols_by_prefix[pref] = {}
            prefix_order.append(pref)
        metric_cols_by_prefix[pref][suf] = c

    if metric_indices:
        first_i = min(metric_indices)
        last_i = max(metric_indices)
        left = cols[:first_i]
        right = cols[last_i + 1 :]
    else:
        # No metric cols found; keep as-is
        left = cols
        right = []
        first_i = None
        last_i = None

    # Build reordered metric block (prefix blocks preserved, suffix order enforced)
    reordered_metrics = []
    for pref in prefix_order:
        d = metric_cols_by_prefix[pref]
        for suf in METRIC_SUFFIX_ORDER:
            if suf in d:
                reordered_metrics.append(d[suf])

    # New header: left + reordered_metrics + right
    if metric_indices:
        new_cols = left + reordered_metrics + right
    else:
        new_cols = cols

    # ---- Fill empties by rules
    # 1) counts -> 0
    for c in df.columns:
        if c.endswith("_REF_COUNT") or c.endswith("_ALT_COUNT"):
            df.loc[df[c].apply(is_empty), c] = "0"

    # 2) numeric CI/MAF -> NA
    for c in df.columns:
        if c.endswith("_MAF") or c.endswith("_LOWER_CI") or c.endswith("_UPPER_CI"):
            df.loc[df[c].apply(is_empty), c] = "NA"

    # 3) IS_GREATER -> "."
    for c in df.columns:
        if c.endswith("_IS_GREATER"):
            df.loc[df[c].apply(is_empty), c] = "."

    # 4) which/platform text flags -> "."
    #    - generic rule: *_which, *_platform
    #    - plus explicitly mentioned ones are covered by this rule anyway
    for c in df.columns:
        if c.endswith("_which") or c.endswith("_platform"):
            df.loc[df[c].apply(is_empty), c] = "."

    # Reorder columns
    df = df[new_cols]

    # Write TSV without quoting
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        out_path,
        sep="\t",
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\",
        na_rep="NA",
    )

    if args.verbose:
        import sys
        sys.stderr.write(f"[INFO] n_rows={df.shape[0]} n_cols={df.shape[1]}\n")
        if first_i is not None:
            sys.stderr.write(f"[INFO] metric_block_header_range={first_i+1}-{last_i+1} (1-based)\n")
        sys.stderr.write(f"[INFO] wrote: {out_path}\n")


if __name__ == "__main__":
    main()


