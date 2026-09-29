#!/usr/bin/env python3
"""
trinuc_denoising.py

Purpose:
    Denoise the variant table with trinucleotide-context-specific ML_QUAL
    thresholds learned in hom_and_single.py. A variant is kept only if its
    ML_QUAL is at or above the threshold of its own triN context AND at or
    above a global minimum cutoff (12.0).

Usage:
    python trinuc_denoising.py \
        --pre <variants.tsv> \
        --thr <trinuc_thresholds.tsv> \
        --out <denoised.tsv>

Inputs:
    --pre   Variant table from conversion.py (must contain 'triN', 'ML_QUAL')
    --thr   Per-context thresholds from hom_and_single.py
            (must contain 'triN', 'best_threshold')

Output:
    --out   Variants passing the filter (same columns as --pre, plus
            'best_threshold' and 'effective_thr')

Notes:
    - Contexts without a learned threshold (e.g. observed in only one
      training class) fall back to the global cutoff.
    - The global cutoff acts as a floor: even if the learned threshold of a
      context is lower than 12.0, variants with ML_QUAL < 12.0 are removed.
    - Variants with missing ML_QUAL or triN never pass (NaN comparisons are
      False, and unmatched triN uses the global cutoff).

Pipeline context:
    Snakemake rule step07c_denoising (inputs: step07a TSV, step07b thresholds).
    The output is converted back to VCF and split into multiread / singleton
    sets in trinucDenoised_integration_refining.py (step07d).
"""

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

    GLOBAL_CUTOFF = 12.0    # minimum ML_QUAL for any variant, and fallback threshold

    sys.stderr.write(f"[INFO] Loading pre-denoise table: {args.pre}\n")
    ffm_df = pd.read_csv(args.pre, sep="\t")

    sys.stderr.write(f"[INFO] Loading trinuc thresholds: {args.thr}\n")
    thr_df = pd.read_csv(args.thr, sep="\t")

    # Sanity check
    for col in ["triN", "ML_QUAL"]:
        if col not in ffm_df.columns:
            sys.exit(f"[ERROR] pre TSV is missing '{col}' column.\n")
    for col in ["triN", "best_threshold"]:
        if col not in thr_df.columns:
            sys.exit(f"[ERROR] thresholds TSV must contain '{col}' column.\n")

    # Attach the context-specific threshold to each variant
    # contexts without a learned threshold use the global cutoff
    merged = ffm_df.merge(thr_df[["triN", "best_threshold"]], on="triN", how="left")
    merged["effective_thr"] = merged["best_threshold"].fillna(GLOBAL_CUTOFF)


    # Keep variants passing BOTH the context-specific threshold and the
    # global floor (the floor matters when a learned threshold is < 12.0)
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
