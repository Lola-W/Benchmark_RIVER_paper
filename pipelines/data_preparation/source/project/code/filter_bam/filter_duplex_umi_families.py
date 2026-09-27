#!/usr/bin/env python3
"""
Filter duplex read families by per-family strand support.

Reads are grouped by UMI family from:
1) `DB:Z:<UMI1>-<UMI2>` (preferred), or
2) read name suffix after the last "_" (fallback), e.g. `<name>_UMI1+UMI2`.

For each family, this script counts:
- `FWD`: reads with FLAG bit 16 unset
- `REV`: reads with FLAG bit 16 set

Families are retained only if:
    FWD >= --min-fwd and REV >= --min-rev

Outputs:
1) Per-family summary table (TSV)
2) Filtered alignment file

The script supports both:
- standard BAM/SAM files via pysam
- headerless SAM-like text partitions (even if file extension is ".bam")
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict
from typing import DefaultDict, Dict, Iterable, Optional, Tuple

try:
    import pysam
except ImportError:  # pragma: no cover
    pysam = None


Counts = DefaultDict[str, list]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count duplex UMI families and remove families below strand thresholds."
    )
    parser.add_argument("input", help="Input alignment file (BAM/SAM or SAM-like text).")
    parser.add_argument("summary_tsv", help="Output TSV with per-UMI strand counts.")
    parser.add_argument(
        "cleaned_output",
        help="Output filtered alignment file. For text inputs, output is text.",
    )
    parser.add_argument(
        "--min-fwd",
        type=int,
        default=1,
        help="Minimum forward reads required per UMI family (default: 1).",
    )
    parser.add_argument(
        "--min-rev",
        type=int,
        default=1,
        help="Minimum reverse reads required per UMI family (default: 1).",
    )
    parser.add_argument(
        "--canonicalize-order",
        action="store_true",
        help="Treat UMI1+UMI2 and UMI2+UMI1 as the same family by sorting UMIs.",
    )
    parser.add_argument(
        "--keep-missing-umi",
        action="store_true",
        help="Keep reads without a parseable UMI in the cleaned output.",
    )
    parser.add_argument(
        "--force-text",
        action="store_true",
        help="Force text parsing mode (for headerless/truncated SAM-like partitions).",
    )
    return parser.parse_args()


def normalize_umi(raw_umi: Optional[str], canonicalize_order: bool) -> Optional[str]:
    if not raw_umi:
        return None
    umi = raw_umi.strip()
    if umi.startswith("DB:Z:"):
        umi = umi[5:]
    umi = umi.replace("-", "+")
    parts = umi.split("+")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        return None
    if canonicalize_order:
        parts = sorted(parts)
    return f"{parts[0]}+{parts[1]}"


def umi_from_qname(qname: str) -> Optional[str]:
    if "_" not in qname:
        return None
    suffix = qname.rsplit("_", 1)[1]
    if "+" in suffix or "-" in suffix:
        return suffix
    return None


def umi_from_sam_fields(qname: str, optional_fields: Iterable[str]) -> Optional[str]:
    for field in optional_fields:
        if field.startswith("DB:Z:"):
            return field[5:]
    return umi_from_qname(qname)


def detect_reader_mode(path: str, force_text: bool) -> str:
    if force_text:
        return "text"
    if pysam is None:
        return "text"

    for mode in ("rb", "r"):
        try:
            with pysam.AlignmentFile(
                path, mode, check_sq=False, require_index=False
            ) as aln:
                it = aln.fetch(until_eof=True)
                try:
                    next(it)
                except StopIteration:
                    pass
            return f"pysam:{mode}"
        except Exception:
            continue
    return "text"


def count_families_text(
    path: str, canonicalize_order: bool
) -> Tuple[Counts, int, int]:
    counts: Counts = defaultdict(lambda: [0, 0])
    total_reads = 0
    missing_umi_reads = 0

    with open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line[0] == "@":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue
            try:
                flag = int(cols[1])
            except ValueError:
                continue

            total_reads += 1
            qname = cols[0]
            raw_umi = umi_from_sam_fields(qname, cols[11:])
            umi = normalize_umi(raw_umi, canonicalize_order)
            if umi is None:
                missing_umi_reads += 1
                continue

            if flag & 16:
                counts[umi][1] += 1
            else:
                counts[umi][0] += 1

    return counts, total_reads, missing_umi_reads


def count_families_pysam(
    path: str, mode: str, canonicalize_order: bool
) -> Tuple[Counts, int, int]:
    counts: Counts = defaultdict(lambda: [0, 0])
    total_reads = 0
    missing_umi_reads = 0

    with pysam.AlignmentFile(path, mode, check_sq=False, require_index=False) as aln:
        for rec in aln.fetch(until_eof=True):
            total_reads += 1
            raw_umi = rec.get_tag("DB") if rec.has_tag("DB") else umi_from_qname(rec.query_name)
            umi = normalize_umi(raw_umi, canonicalize_order)
            if umi is None:
                missing_umi_reads += 1
                continue

            if rec.flag & 16:
                counts[umi][1] += 1
            else:
                counts[umi][0] += 1

    return counts, total_reads, missing_umi_reads


def write_summary(
    summary_tsv: str, counts: Dict[str, list], min_fwd: int, min_rev: int
) -> None:
    out_dir = os.path.dirname(summary_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    def sort_key(item: Tuple[str, list]) -> Tuple[int, int, str]:
        umi, (fwd, rev) = item
        return (-(fwd + rev), -fwd, umi)

    with open(summary_tsv, "wt", encoding="utf-8") as out:
        out.write("UMI_UMI\tFWD\tREV\tFWD+REV\tPASS\n")
        for umi, (fwd, rev) in sorted(counts.items(), key=sort_key):
            keep = int(fwd >= min_fwd and rev >= min_rev)
            out.write(f"{umi}\t{fwd}\t{rev}\t{fwd}+{rev}\t{keep}\n")


def filter_text(
    input_path: str,
    output_path: str,
    passing_families: set,
    canonicalize_order: bool,
    keep_missing_umi: bool,
) -> Tuple[int, int, int]:
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    kept = 0
    removed = 0
    missing_umi_reads = 0

    with open(input_path, "rt", encoding="utf-8", errors="replace") as src, open(
        output_path, "wt", encoding="utf-8"
    ) as dst:
        for line in src:
            if not line or line[0] == "@":
                dst.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                dst.write(line)
                continue

            raw_umi = umi_from_sam_fields(cols[0], cols[11:])
            umi = normalize_umi(raw_umi, canonicalize_order)
            if umi is None:
                missing_umi_reads += 1
                if keep_missing_umi:
                    dst.write(line)
                    kept += 1
                else:
                    removed += 1
                continue

            if umi in passing_families:
                dst.write(line)
                kept += 1
            else:
                removed += 1

    return kept, removed, missing_umi_reads


def pysam_output_mode(cleaned_output: str) -> str:
    lower = cleaned_output.lower()
    if lower.endswith(".sam"):
        return "w"
    return "wb"


def filter_pysam(
    input_path: str,
    output_path: str,
    read_mode: str,
    passing_families: set,
    canonicalize_order: bool,
    keep_missing_umi: bool,
) -> Tuple[int, int, int]:
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    kept = 0
    removed = 0
    missing_umi_reads = 0

    with pysam.AlignmentFile(
        input_path, read_mode, check_sq=False, require_index=False
    ) as src:
        out_mode = pysam_output_mode(output_path)
        with pysam.AlignmentFile(output_path, out_mode, header=src.header) as dst:
            for rec in src.fetch(until_eof=True):
                raw_umi = rec.get_tag("DB") if rec.has_tag("DB") else umi_from_qname(rec.query_name)
                umi = normalize_umi(raw_umi, canonicalize_order)
                if umi is None:
                    missing_umi_reads += 1
                    if keep_missing_umi:
                        dst.write(rec)
                        kept += 1
                    else:
                        removed += 1
                    continue

                if umi in passing_families:
                    dst.write(rec)
                    kept += 1
                else:
                    removed += 1

    return kept, removed, missing_umi_reads


def main() -> int:
    args = parse_args()
    if args.min_fwd < 0 or args.min_rev < 0:
        print("ERROR: --min-fwd and --min-rev must be >= 0.", file=sys.stderr)
        return 2

    reader_mode = detect_reader_mode(args.input, args.force_text)
    print(f"[INFO] Reader mode: {reader_mode}", file=sys.stderr)

    if reader_mode.startswith("pysam:"):
        read_mode = reader_mode.split(":", 1)[1]
        counts, total_reads, missing_umi_count = count_families_pysam(
            args.input, read_mode, args.canonicalize_order
        )
    else:
        counts, total_reads, missing_umi_count = count_families_text(
            args.input, args.canonicalize_order
        )

    passing_families = {
        umi
        for umi, (fwd, rev) in counts.items()
        if fwd >= args.min_fwd and rev >= args.min_rev
    }

    write_summary(args.summary_tsv, counts, args.min_fwd, args.min_rev)

    if reader_mode.startswith("pysam:"):
        read_mode = reader_mode.split(":", 1)[1]
        kept, removed, missing_umi_filter = filter_pysam(
            args.input,
            args.cleaned_output,
            read_mode,
            passing_families,
            args.canonicalize_order,
            args.keep_missing_umi,
        )
    else:
        kept, removed, missing_umi_filter = filter_text(
            args.input,
            args.cleaned_output,
            passing_families,
            args.canonicalize_order,
            args.keep_missing_umi,
        )

    print(
        "[INFO] Total reads scanned: "
        f"{total_reads}; families found: {len(counts)}; "
        f"families passing {args.min_fwd}+{args.min_rev}: {len(passing_families)}",
        file=sys.stderr,
    )
    print(
        "[INFO] Missing-UMI reads during counting: "
        f"{missing_umi_count}; during filtering: {missing_umi_filter}",
        file=sys.stderr,
    )
    print(
        f"[INFO] Cleaned output written: {args.cleaned_output} (kept={kept}, removed={removed})",
        file=sys.stderr,
    )
    print(f"[INFO] Summary written: {args.summary_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
