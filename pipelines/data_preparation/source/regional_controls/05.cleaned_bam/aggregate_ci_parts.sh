#!/usr/bin/env bash
set -euo pipefail

IN_DIR="$1"
OUT_TSV="$2"

first="$(ls -1 "${IN_DIR}"/part_*.tsv | sort | head -n 1)"

mkdir -p "$(dirname "$OUT_TSV")"

{
  head -n 1 "$first"
  for f in $(ls -1 "${IN_DIR}"/part_*.tsv | sort); do
    tail -n +2 "$f"
  done
} > "$OUT_TSV"

