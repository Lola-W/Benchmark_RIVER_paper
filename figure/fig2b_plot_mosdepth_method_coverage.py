#!/usr/bin/env python3
"""Plot mosdepth method-coverage curves with the shared NBT figure style.

The data reader and axis geometry follow the original reproducibility script.
This version standardizes method names and colors, uses Arial, omits the
automatic header, and writes a transparent vector PDF by default.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from statistics import fmean, median
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        # Keep text editable and searchable in PDF/PS output.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


DEFAULT_CDF_SUFFIX = ".keep.global.cdf.tsv"
DEFAULT_PDF_SUFFIX = ".keep.global.pdf.tsv"

SUBMISSION_ROOT = Path(__file__).resolve().parent / ".."

# Mirrors R/nature_biotech_ggplot.R so a method has the same identity in every
# figure. Main curves use the exact stable color; raw curves use a lighter shade.
NBT_METHOD_COLORS = {
    "ppmSeq": "#0072B2",
    "UDSeq": "#009E73",
    "NanoSeq": "#E69F00",
    "HiDEF-seq": "#D55E00",
    "Illumina": "#CC79A7",
    "PTA": "#56B4E9",
}

METHOD_NAMES = {
    "hidefseq": "HiDEF-seq",
    "illumina": "Illumina",
    "nanoseq": "NanoSeq",
    "ppmseq": "ppmSeq",
    "pta": "PTA",
    "udseq": "UDSeq",
}


def default_input_dir() -> Path:
    """Use ../raw_data/coverage, with the existing environment override."""
    local = SUBMISSION_ROOT / "raw_data" / "coverage"
    if local.is_dir():
        return local

    configured_root = os.environ.get("MOSDEPTH_REPRO_ROOT")
    if configured_root:
        configured = Path(configured_root).expanduser() / "data" / "coverage"
        if configured.is_dir():
            return configured

    return local


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot mosdepth coverage curves in the shared NBT style."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir(),
        help="Coverage directory containing method folders "
        f"(default: {default_input_dir()})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SUBMISSION_ROOT / "output" / "method_coverage_cdf_0_100_nbt.pdf",
        help="Output path; transparent vector PDF is the default.",
    )
    parser.add_argument(
        "--kind",
        choices=("cdf", "pdf"),
        default="cdf",
        help="Use `.keep.global.cdf.tsv` (default) or `.keep.global.pdf.tsv`.",
    )
    parser.add_argument(
        "--agg",
        choices=("mean", "median"),
        default="mean",
        help="Statistic used when a method folder has multiple files.",
    )
    parser.add_argument("--xmin", type=int, default=0)
    parser.add_argument("--xmax", type=int, default=100)
    parser.add_argument(
        "--include-tmp",
        action="store_true",
        help="Include the `tmp` data folder (excluded by default).",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional title. No title/header is drawn by default.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution used only for raster output.",
    )
    parser.add_argument(
        "--opaque",
        action="store_true",
        help="Use an opaque white background instead of transparency.",
    )
    return parser.parse_args()


def expected_suffix(kind: str) -> str:
    if kind == "cdf":
        return DEFAULT_CDF_SUFFIX
    if kind == "pdf":
        return DEFAULT_PDF_SUFFIX
    raise ValueError(f"Unsupported kind: {kind}")


def discover_method_files(
    input_dir: Path, suffix: str, include_tmp: bool
) -> dict[str, list[Path]]:
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    methods: dict[str, list[Path]] = {}
    for child in sorted(input_dir.iterdir(), key=lambda path: path.name.lower()):
        if not child.is_dir() or child.name.startswith("."):
            continue
        if child.name == "tmp" and not include_tmp:
            continue

        files = sorted(child.glob(f"*{suffix}"))
        if files:
            methods[child.name] = files
        else:
            print(f"[skip] {child.name}: no files matching *{suffix}", file=sys.stderr)
    return methods


def read_two_col_tsv(path: Path) -> dict[int, float]:
    curve: dict[int, float] = {}
    with path.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(f"{path}: line {line_number} has fewer than 2 columns")
            try:
                coverage = int(row[0])
                proportion = float(row[1])
            except ValueError as exc:
                raise ValueError(
                    f"{path}: line {line_number} is not numeric: {row[:2]}"
                ) from exc
            curve[coverage] = proportion
    if not curve:
        raise ValueError(f"{path}: no data rows found")
    return curve


def aggregate_values(values: Iterable[float], agg: str) -> float:
    vals = list(values)
    if not vals:
        raise ValueError("No values to aggregate")
    if agg == "mean":
        return fmean(vals)
    if agg == "median":
        return median(vals)
    raise ValueError(f"Unsupported agg: {agg}")


def aggregate_method_curves(
    files: list[Path], x_values: list[int], agg: str
) -> list[float]:
    curves = [read_two_col_tsv(path) for path in files]
    values: list[float] = []
    for coverage in x_values:
        value = aggregate_values(
            (curve.get(coverage, 0.0) for curve in curves), agg=agg
        )
        values.append(min(max(value, 0.0), 1.0))
    return values


def split_method_variant(method: str) -> tuple[str, str]:
    method_key = method.lower()
    if method_key.endswith("_raw"):
        return method_key[: -len("_raw")], "raw"
    return method_key, "main"


def display_family(method_key: str) -> str:
    return METHOD_NAMES.get(method_key, method_key)


def display_label(method: str, file_count: int, agg: str) -> str:
    family_key, variant = split_method_variant(method)
    label = display_family(family_key)
    if variant == "raw":
        label = f"{label} (raw)"
    if file_count > 1:
        label = f"{label} ({agg}; n = {file_count})"
    return label


def lighten(color: str, amount: float = 0.42) -> tuple[float, float, float]:
    rgb = to_rgb(color)
    return tuple(channel + (1.0 - channel) * amount for channel in rgb)


def method_appearance(method: str) -> tuple[str | tuple[float, float, float], str]:
    family_key, variant = split_method_variant(method)
    family = display_family(family_key)
    base_color = NBT_METHOD_COLORS.get(family, "#4D4D4D")
    if variant == "raw":
        return lighten(base_color), "--"
    return base_color, "-"


def main() -> int:
    args = parse_args()
    if args.xmin > args.xmax:
        raise SystemExit("--xmin must be <= --xmax")

    suffix = expected_suffix(args.kind)
    methods = discover_method_files(
        args.input_dir, suffix=suffix, include_tmp=args.include_tmp
    )
    if not methods:
        raise SystemExit(
            f"No method folders with files matching *{suffix} under {args.input_dir}"
        )

    x_values = list(range(args.xmin, args.xmax + 1))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    transparent = not args.opaque
    if transparent:
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    for method, files in methods.items():
        y_values = aggregate_method_curves(files, x_values=x_values, agg=args.agg)
        color, line_style = method_appearance(method)
        ax.plot(
            x_values,
            y_values,
            linewidth=2,
            linestyle=line_style,
            color=color,
            label=display_label(method, len(files), args.agg),
        )

    if args.kind == "cdf":
        ylabel = "Genome fraction at or above depth"
        ax.set_ylim(-0.02, 1.02)
    else:
        ylabel = "Genome fraction at depth"
        ax.set_ylim(bottom=-0.02)

    if args.title:
        ax.set_title(args.title)
    ax.set_xlabel("Coverage depth (×)")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3, linewidth=0.8)
    ax.legend(frameon=False)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, transparent=transparent)
    plt.close(fig)

    print(f"Saved figure: {args.output}")
    print(f"Input dir: {args.input_dir}")
    print(f"Methods plotted: {', '.join(methods)}")
    for method, files in methods.items():
        if len(files) > 1:
            print(f"  {method}: aggregated {len(files)} files using {args.agg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
