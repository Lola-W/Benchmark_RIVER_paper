#!/usr/bin/env python3
"""
hom_and_single.py

Purpose:
    Learn a trinucleotide-context-specific ML_QUAL threshold. 
    Two labeled training sets are combined:
        - FP-like (label 0): outMap variants shared with sbsMap, i.e. loci supported by a single read (singletons, likely artifacts)
        - TP-like (label 1): outMap variants shared with homMap, i.e. loci where most reads support the variant (likely homozygous germline)
    For each trinucleotide context, the ML_QUAL cutoff that best separates the two classes (Youden's J) is written to a table.
    The approach follows the trinucleotide motif-specific thresholding.

Usage:
    python hom_and_single.py \
        --single <outMap.intersect_sbsMap.vcf.gz> \
        --homo   <outMap.intersect_homMap.vcf.gz> \
        --outdir <output_dir> \
        --thresh-out <trinuc_thresholds.tsv>

Inputs:
    --single      FP-like VCF (output of overlap_sbs.sh)
    --homo        TP-like VCF (output of overlap_hom.sh)

Outputs:
    --thresh-out  Per-context table: triN, AUC, best_threshold,
                  sensitivity, specificity
    --outdir      Also receives training_dataset.tsv (the labeled training
                  data actually used: ML_QUAL, triN, label)

Notes:
    - A VCF with >= 10,000,000 records is randomly subsampled to exactly 1% (seed 42) to limit memory use.
    - triN is the trinucleotide substitution context in the orientation in which the read was sequenced (192 motifs).
    - ML_QUAL is the raw (pre-recalibration) score of the SRSNV classifier, stored in INFO.
      Thresholds are on the ML_QUAL scale. Higher values are assumed to indicate a more likely true variant.
      
Pipeline context:
    Snakemake rule step07b_training. The thresholds table is consumed by
    trinuc_denoising.py (step07c).
"""

import gzip
import os
import random
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from sklearn.metrics import roc_curve, roc_auc_score


def rev_triN(triN):
    """Reverse complement of a trinucleotide context.

    Format is "X[Y>Z]W" (e.g. "A[C>T]G" -> "C[G>A]T"). Same as in conversion.py.
    """
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return comp[triN[6]] + '[' + comp[triN[2]] + '>' + comp[triN[4]] + ']' + comp[triN[0]]


def count_records(filename):
    """Count non-header records in a (gzipped) VCF without parsing them."""
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as f:
        return sum(1 for line in f if not line.startswith("#"))


def vcf_to_df(filename, label, frac=0.01, threshold=10_000_000):
    """Read a ppmSeq featuremap VCF into a DataFrame.

    VCFs with >= `threshold` records are subsampled to exactly `frac`
    (default 1%) to limit memory; smaller VCFs are used in full.
    Returns columns ML_QUAL, triN and a constant class label
    (0 = FP-like, 1 = TP-like).
    """
    opener = gzip.open if filename.endswith(".gz") else open

    # Pass 1: count records, then choose which record indices to keep.
    # Only the chosen indices are stored, not the lines themselves.
    n_total = count_records(filename)
    keep = None
    if n_total >= threshold:
        keep = set(random.sample(range(n_total), int(n_total * frac)))
    sys.stderr.write(
        f"[INFO] {filename}: {n_total} records, "
        f"{'subsampled to ' + str(len(keep)) if keep is not None else 'using all'}\n"
    )

    ml_quals = []
    triNs = []

    # Pass 2: parse only the selected records
    with opener(filename, "rt") as f:
        i = -1
        for line in f:
            if line.startswith("#"):
                continue
            i += 1
            if keep is not None and i not in keep:
                continue

            fields = line.strip().split("\t")  # CHROM POS ID REF ALT QUAL FILTER INFO
            info_dict = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)

            # Strand-normalized trinucleotide context (see conversion.py)
            v = info_dict.get("trinuc_context_with_alt", "")
            triN = f"{v[0]}[{v[1]}>{v[3]}]{v[2]}" if len(v) >= 4 else None
            if info_dict.get("X_FLAGS") in ("16", "1040") and triN:
                triN = rev_triN(triN)  # reverse-strand read

            ml_quals.append(info_dict.get("ML_QUAL"))
            triNs.append(triN)

    df = pd.DataFrame({"ML_QUAL": ml_quals, "triN": triNs, "label": label})
    df["ML_QUAL"] = pd.to_numeric(df["ML_QUAL"], errors='coerce').astype("float32")
    return df


def compute_best_threshold(y_true, y_score):
    """Find the ML_QUAL threshold maximizing Youden's J (TPR - FPR).

    y_true: class labels (0 or 1); y_score: ML_QUAL.
    Returns (threshold, sensitivity, specificity) at the optimum.
    """
    fpr, tpr, thr = roc_curve(y_true, y_score)
    idx = np.argmax(tpr - fpr)
    return thr[idx], tpr[idx], 1 - fpr[idx]


def main():
    random.seed(42)  # fixed seed so the subsample is reproducible
    parser = ArgumentParser(
        prog='hom_and_single.py',
        description='Combine FP-like (sbsMap) + TP-like (homMap) VCFs to compute per-trinuc ML_QUAL thresholds.'
    )
    parser.add_argument("--single", required=True, help="FP-like intersect_sbsMap .vcf.gz")
    parser.add_argument("--homo",   required=True, help="TP-like intersect_homMap .vcf.gz")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--thresh-out", required=True, help="Output trinuc thresholds TSV path")
    parser.add_argument("--subsample-frac", type=float, default=0.01,
                        help="Fraction kept for large VCFs (default: 0.01)")
    parser.add_argument("--subsample-threshold", type=int, default=10_000_000,
                        help="VCFs with >= this many records are subsampled (default: 10,000,000)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    sys.stderr.write(f"[INFO] Loading FP-like VCF: {args.single}\n")
    single_df = vcf_to_df(args.single, label=0, frac=args.subsample_frac, threshold=args.subsample_threshold)

    sys.stderr.write(f"[INFO] Loading TP-like VCF: {args.homo}\n")
    homo_df = vcf_to_df(args.homo, label=1, frac=args.subsample_frac, threshold=args.subsample_threshold)

    # Combine both classes and save the training dataset
    combined = pd.concat([single_df, homo_df], ignore_index=True)
    training_tsv = os.path.join(args.outdir, "training_dataset.tsv")
    combined.to_csv(training_tsv, sep="\t", index=False)
    sys.stderr.write(f"[INFO] Wrote: {training_tsv}\n")

    # Compute one threshold per trinucleotide context
    records = []
    for triN, sub in combined.groupby("triN"):
        if sub["label"].nunique() < 2:
            continue  # both classes are needed to compute ROC
        ml_qual = sub["ML_QUAL"].dropna()
        labels  = sub["label"].loc[ml_qual.index]
        auc = roc_auc_score(labels, ml_qual)
        best_thr, sens, spec = compute_best_threshold(labels, ml_qual)
        records.append((triN, auc, best_thr, sens, spec))

    pd.DataFrame(records, columns=["triN","AUC","best_threshold","sensitivity","specificity"])\
      .to_csv(args.thresh_out, sep="\t", index=False)
    sys.stderr.write(f"[INFO] Wrote: {args.thresh_out}\n")


if __name__ == "__main__":
    main()
