#!/usr/bin/env python3

import gzip
import os
import random
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from sklearn.metrics import roc_curve, roc_auc_score


def rev_triN(triN):
    "Reverse complement of a trinucleotide context."
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return comp[triN[6]] + '[' + comp[triN[2]] + '>' + comp[triN[4]] + ']' + comp[triN[0]]


def vcf_to_df(filename, label):
    """Read a ppmSeq featuremap VCF into a DataFrame with ML_QUAL, triN, and label."""
    opener = gzip.open if filename.endswith(".gz") else open

    ml_quals = []
    triNs = []

    with opener(filename, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if random.random() > 0.01: #1%
                continue
            fields = line.strip().split("\t")
            info_dict = dict(item.split("=", 1) for item in fields[7].split(";") if "=" in item)

            v = info_dict.get("trinuc_context_with_alt", "")
            triN = f"{v[0]}[{v[1]}>{v[3]}]{v[2]}" if len(v) >= 4 else None
            if info_dict.get("X_FLAGS") in ("16", "1040") and triN:
                triN = rev_triN(triN)

            ml_quals.append(info_dict.get("ML_QUAL"))
            triNs.append(triN)

    df = pd.DataFrame({"ML_QUAL": ml_quals, "triN": triNs, "label": label})
    df["ML_QUAL"] = pd.to_numeric(df["ML_QUAL"], errors='coerce').astype("float32")
    return df


def compute_best_threshold(y_true, y_score): # y_ture: label(0 or 1), y_score: ML_QUAL
    """Find best ML_QUAL threshold via Youden's J (max TPR - FPR)."""
    fpr, tpr, thr = roc_curve(y_true, y_score)
    idx = np.argmax(tpr - fpr)
    return thr[idx], tpr[idx], 1 - fpr[idx]


def main():
    random.seed(42)  ## Fix seed
    parser = ArgumentParser(
        prog='hom_and_single.py',
        description='Combine FP-like (sbsMap) + TP-like (homMap) VCFs to compute per-trinuc ML_QUAL thresholds.'
    )
    parser.add_argument("--single", required=True, help="FP-like intersect_sbsMap .vcf.gz")
    parser.add_argument("--homo",   required=True, help="TP-like intersect_homMap .vcf.gz")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--thresh-out", required=True, help="Output trinuc thresholds TSV path")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    sys.stderr.write(f"[INFO] Loading FP-like VCF: {args.single}\n")
    single_df = vcf_to_df(args.single, label=0)  # FP-like: label 0

    sys.stderr.write(f"[INFO] Loading TP-like VCF: {args.homo}\n")
    homo_df = vcf_to_df(args.homo, label=1)       # TP-like: label 1

    # Combine and write training dataset
    combined = pd.concat([single_df, homo_df], ignore_index=True)
    training_tsv = os.path.join(args.outdir, "training_dataset.tsv")
    combined.to_csv(training_tsv, sep="\t", index=False)
    sys.stderr.write(f"[INFO] Wrote: {training_tsv}\n")

    # Compute per-trinuc ML_QUAL thresholds
    records = []
    for triN, sub in combined.groupby("triN"):
        if sub["label"].nunique() < 2:
            continue  # Skip if only one class present
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
