#!/usr/bin/env python3
"""Assemble the final plotting matrix from historical 32-context count files.

New submission utility, not the historical merger.sh. Counts are joined by
context and columns explicitly mapped to the final plotting names.
Python 3.9+, standard library only. Run on a compute node on TSCC.
"""
import argparse
import csv
from pathlib import Path

CONTEXTS = tuple(a + b + c for b in "CT" for a in "ACGT" for c in "ACGT")
DEFAULT_MAPPING = Path(__file__).resolve().parents[1] / "config/trinucleotide_columns.tsv"

def read_mapping(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["plot_column", "counts_file"]:
            raise ValueError("Mapping header must be plot_column<TAB>counts_file")
        mapping = list(reader)
    if not mapping:
        raise ValueError("Empty column mapping")
    for row in mapping:
        if set(row) != {"plot_column", "counts_file"} or not all(row.values()):
            raise ValueError("Malformed column mapping row")
        filename = Path(row["counts_file"])
        if filename.name != row["counts_file"] or filename.name in {".", ".."}:
            raise ValueError("counts_file must be a filename within --counts-dir")
    names = [row["plot_column"] for row in mapping]
    if len(set(names)) != len(names) or "trinuc" in names:
        raise ValueError("Duplicate or reserved plotting column")
    return mapping

def read_counts(path):
    counts = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) != 2 or header[0] != "trinuc":
            raise ValueError(f"Expected a two-column trinuc/count table: {path}")
        for line, row in enumerate(reader, 2):
            if len(row) != 2:
                raise ValueError(f"Expected two fields: {path}:{line}")
            context, value = row
            if context in counts or context not in CONTEXTS:
                raise ValueError(f"Duplicate or invalid context {context!r}: {path}:{line}")
            try:
                count = int(value)
            except ValueError as exc:
                raise ValueError(f"Noninteger count: {path}:{line}") from exc
            if count < 0:
                raise ValueError(f"Negative count: {path}:{line}")
            counts[context] = count
    if set(counts) != set(CONTEXTS):
        raise ValueError(f"Expected all 32 pyrimidine contexts: {path}")
    return counts

def assemble(counts_dir, mapping_path, output):
    mapping = read_mapping(mapping_path)
    tables = [read_counts(counts_dir / row["counts_file"]) for row in mapping]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trinuc"] + [row["plot_column"] for row in mapping])
        for context in CONTEXTS:
            writer.writerow([context] + [table[context] for table in tables])
    return len(tables)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        columns = assemble(args.counts_dir, args.mapping, args.output)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{exc}\n")
    print(f"Wrote 32 contexts x {columns} profiles to {args.output}")

if __name__ == "__main__":
    main()
