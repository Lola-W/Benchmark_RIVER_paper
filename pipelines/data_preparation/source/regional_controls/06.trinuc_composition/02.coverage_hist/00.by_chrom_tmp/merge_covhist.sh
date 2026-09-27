#!/usr/bin/env bash
set -euo pipefail

OUTDIR="../01.merged"
mkdir -p "$OUTDIR"

platforms=$(
  ls *_coverage_hist.chr*.tsv \
  | sed -E 's/\.chr[0-9XY]+\.tsv$//' \
  | sort -u
)

for p in $platforms; do
  out="${OUTDIR}/${p}.merged.tsv"
  tmp="$(mktemp)"

  # 1) concat chr1-22,X
  # 2) drop per-file headers robustly
  # 3) sum by depth
  # 4) output ONLY numeric lines to tmp
  cat ${p}.chr{1..22}.tsv ${p}.chrX.tsv 2>/dev/null \
  | awk '
      ($1=="depth" && $2=="count"){next}
      ($1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/){cnt[$1]+=$2}
      END{for (d in cnt) print d "\t" cnt[d]}
    ' > "$tmp"

  # sort numeric ONLY, then prepend header
  {
    printf "depth\tcount\n"
    sort -k1,1n "$tmp"
  } > "$out"

  rm -f "$tmp"
  echo "[OK] $out"
done
