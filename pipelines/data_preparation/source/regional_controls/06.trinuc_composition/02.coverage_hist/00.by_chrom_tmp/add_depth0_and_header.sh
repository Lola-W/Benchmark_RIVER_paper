#!/usr/bin/env bash
set -euo pipefail

REF="ref_chr_count.tsv"
mkdir -p tmp

# Load ref counts into associative array
declare -A REFCOUNT
while IFS=$'\t' read -r chr cnt; do
  [[ -z "${chr}" ]] && continue
  REFCOUNT["$chr"]="$cnt"
done < "$REF"

# Process all coverage hist files
shopt -s nullglob
for f in *_coverage_hist.chr*.tsv; do
  base="$(basename "$f")"

  # Extract chromosome from filename: ... .chr1.tsv / .chrX.tsv
  chr="$(sed -n 's/.*\.\(chr[0-9XY]\+\)\.tsv/\1/p' <<< "$base")"
  if [[ -z "${chr}" ]]; then
    echo "[SKIP] cannot parse chr from: $f" >&2
    continue
  fi

  ref="${REFCOUNT[$chr]:-}"
  if [[ -z "${ref}" ]]; then
    echo "[SKIP] no ref count for $chr (file: $f)" >&2
    continue
  fi

  # Sum 2nd column
  sum="$(awk '{s+=$2} END{printf "%.0f\n", s+0}' "$f")"

  # depth0 = ref - sum
  depth0="$(( ref - sum ))"
  if (( depth0 < 0 )); then
    echo "[WARN] depth0 < 0 for $f (ref=$ref sum=$sum depth0=$depth0)" >&2
  fi

  # Move original to tmp/
  mv -f "$f" "tmp/$base"

  # Write new file with header + depth0 row + original content appended
  {
    printf "depth\tcount\n"
    printf "0\t%d\n" "$depth0"
    cat "tmp/$base"
  } > "$f"

  echo -e "[OK]\t$base\tchr=$chr\tref=$ref\tsum=$sum\tdepth0=$depth0"
done
