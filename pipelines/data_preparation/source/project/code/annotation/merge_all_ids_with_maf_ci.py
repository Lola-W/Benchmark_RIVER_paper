#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd


def chrom_to_variant_prefix(chrom: str) -> str:
    # chr1 -> 1, chrX -> X
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def load_maf_ci(maf_path: Path, prefix: str) -> pd.DataFrame:
    """
    Load one maf_ci.tsv and return a DF indexed by variant_id with prefixed metric columns.
    """
    df = pd.read_csv(maf_path, sep="\t", dtype=str)

    # Some files have '#CHROM', some might have 'CHROM'
    chrom_col = "#CHROM" if "#CHROM" in df.columns else ("CHROM" if "CHROM" in df.columns else None)
    if chrom_col is None:
        raise ValueError(f"Cannot find CHROM column in {maf_path}. Columns={list(df.columns)[:20]}")

    required = [chrom_col, "POS", "REF", "ALT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {maf_path}")

    # Build variant_id = chrom(without chr)-POS-REF-ALT
    df["variant_id"] = (
        df[chrom_col].map(chrom_to_variant_prefix)
        + "-"
        + df["POS"]
        + "-"
        + df["REF"]
        + "-"
        + df["ALT"]
    )

    # Drop location cols (already represented by variant_id)
    drop_cols = [chrom_col, "POS", "REF", "ALT"]
    metric_cols = [c for c in df.columns if c not in drop_cols + ["variant_id"]]

    # Keep only variant_id + metrics
    out = df[["variant_id"] + metric_cols].copy()

    # Prefix metric columns
    out = out.rename(columns={c: f"{prefix}_{c}" for c in metric_cols})

    # Deduplicate if needed (keep first)
    out = out.drop_duplicates(subset=["variant_id"], keep="first")

    # Index by variant_id for fast joins
    out = out.set_index("variant_id")

    return out


def build_maf_table(bam_pileup_dir: Path) -> pd.DataFrame:
    """
    Load all platforms' maf_ci.tsv and merge into one big metrics table keyed by variant_id.
    """
    tables = []

    # One-file platforms
    platform_files = {
        "illumina": bam_pileup_dir / "illumina" / "baseQ0.union.all.snv_only.maf_ci.tsv",
        "udseq": bam_pileup_dir / "udseq" / "baseQ0.union.all.snv_only.maf_ci.tsv",
        "nanoseq": bam_pileup_dir / "nanoseq" / "baseQ0.union.all.snv_only.maf_ci.tsv",
        "hidefseq": bam_pileup_dir / "hidefseq" / "baseQ0.union.all.snv_only.maf_ci.tsv",
        "ppmseq": bam_pileup_dir / "ppmseq" / "baseQ0.union.all.snv_only.maf_ci.tsv",
    }

    for platform, f in platform_files.items():
        if not f.exists():
            raise FileNotFoundError(f"Missing expected file: {f}")
        tables.append(load_maf_ci(f, platform))

    # Single-cell: many files, make per-neuron prefixed columns
    sc_dir = bam_pileup_dir / "single_cell"
    sc_files = sorted(sc_dir.glob("*.maf_ci.tsv"))
    if not sc_files:
        raise FileNotFoundError(f"No single-cell maf_ci files found under: {sc_dir}")

    for f in sc_files:
        # Example: 7614_single_neuron_10.maf_ci.tsv -> single_cell_7614_single_neuron_10
        sample = f.name.replace(".maf_ci.tsv", "")
        prefix = f"single_cell_{sample}"
        tables.append(load_maf_ci(f, prefix))

    # Merge all metric tables on variant_id (outer so we keep any variant seen in any platform)
    merged = tables[0]
    for t in tables[1:]:
        merged = merged.join(t, how="outer")

    return merged


def merge_in_chunks(all_ids_tsv: Path, maf_table: pd.DataFrame, out_tsv: Path, chunksize: int):
    """
    Stream through all_ids_df in chunks and left-join maf_table by variant_id, write one huge TSV.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    first = True
    reader = pd.read_csv(all_ids_tsv, sep="\t", dtype=str, chunksize=chunksize)

    for chunk in reader:
        if "variant_id" not in chunk.columns:
            raise ValueError(f"'variant_id' not found in {all_ids_tsv}")

        chunk = chunk.set_index("variant_id")
        chunk = chunk.join(maf_table, how="left")
        chunk = chunk.reset_index()

        chunk.to_csv(out_tsv, sep="\t", index=False, mode="w" if first else "a", header=first)
        first = False


def main():
    ap = argparse.ArgumentParser(
        description="Merge all_ids_df.annotated.full.tsv with per-platform maf_ci.tsv files (prefixed columns)."
    )
    ap.add_argument("--ids", required=True, type=Path,
                    help="Path to all_ids_df.annotated.full.tsv")
    ap.add_argument("--bamdir", required=True, type=Path,
                    help="Path to bam_pileup/result directory (contains illumina/, udseq/, etc)")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output merged TSV")
    ap.add_argument("--chunksize", type=int, default=200_000,
                    help="Rows per chunk when streaming all_ids_df (default: 200k)")
    args = ap.parse_args()

    maf_table = build_maf_table(args.bamdir)
    merge_in_chunks(args.ids, maf_table, args.out, args.chunksize)
    print(f"[OK] Wrote merged TSV: {args.out}")


if __name__ == "__main__":
    main()


# python ./work/project/code/annotation/merge_all_ids_with_maf_ci.py \
#   --ids ./work/project/data/annotation/all_ids_df.annotated.full.tsv \
#   --bamdir ./work/project/data/bam_pileup/result \
#   --out ./work/project/data/annotation/all_ids_df.annotated.full.with_maf_ci.tsv