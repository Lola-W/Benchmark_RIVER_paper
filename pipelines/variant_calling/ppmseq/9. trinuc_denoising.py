#!/usr/bin/env python3
import sys
import pandas as pd
from argparse import ArgumentParser


def main():
    parser = ArgumentParser(
        prog="trinuc_denoising.py",
        description="Denoise ppmSeq reads using triN-specific ML_QUAL thresholds, "
                    "with a global minimum cutoff enforced regardless of learned threshold."
    )
    parser.add_argument("--pre", required=True, help="Pre-denoise TSV. Must contain 'triN' and 'ML_QUAL'.")
    parser.add_argument("--thr", required=True, help="trinuc_thresholds.tsv. Must contain 'triN' and 'best_threshold'.")
    parser.add_argument("--out", required=True, help="Output denoised TSV path.")
    args = parser.parse_args()

    GLOBAL_CUTOFF = 12.0

    sys.stderr.write(f"[INFO] Loading pre-denoise table: {args.pre}\n")
    ffm_df = pd.read_csv(args.pre, sep="\t")

    sys.stderr.write(f"[INFO] Loading trinuc thresholds: {args.thr}\n")
    thr_df = pd.read_csv(args.thr, sep="\t")

    # sanity check (v2에서 빠져있던 컬럼 검증도 같이 복원)
    for col in ["triN", "ML_QUAL"]:
        if col not in ffm_df.columns:
            sys.exit(f"[ERROR] pre TSV is missing '{col}' column.\n")
    for col in ["triN", "best_threshold"]:
        if col not in thr_df.columns:
            sys.exit(f"[ERROR] thresholds TSV must contain '{col}' column.\n")

    # Merge triN-specific thresholds; use GLOBAL_CUTOFF where threshold is missing
    merged = ffm_df.merge(thr_df[["triN", "best_threshold"]], on="triN", how="left")
    merged["effective_thr"] = merged["best_threshold"].fillna(GLOBAL_CUTOFF)

    # Keep variants where ML_QUAL passes BOTH:
    #   - its triN-specific (or global-fallback) threshold
    #   - the GLOBAL_CUTOFF floor, even if best_threshold was learned lower than this
    # (v2는 이 GLOBAL_CUTOFF 하한선 조건이 누락되어 있었음 -> best_threshold가 12.0보다
    #  낮게 학습된 triN context에서 통과 레코드가 11.57배까지 급증하는 원인이었음)
    keep = (
        (merged["ML_QUAL"] >= merged["effective_thr"]) &
        (merged["ML_QUAL"] >= GLOBAL_CUTOFF)
    )

    denoised = merged[keep].reset_index(drop=True)
    denoised.to_csv(args.out, sep="\t", index=False)

    n_total = len(merged)
    n_kept = len(denoised)
    sys.stderr.write(f"[INFO] Total input records: {n_total}\n")
    sys.stderr.write(f"[INFO] Kept records:        {n_kept}\n")
    sys.stderr.write(f"[INFO] Wrote: {args.out}\n")


if __name__ == "__main__":
    main()
