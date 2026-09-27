#!/usr/bin/env bash
set -euo pipefail

INDIR="./work/project/data/mutational_signatures/sigprofiler/caller_split"
OUTDIR="./work/project/data/mutational_signatures/sigprofiler/caller_split_noputative_final"

mkdir -p "$OUTDIR"

for f in "$INDIR"/*.tsv; do
  base="$(basename "$f")"
  out="$OUTDIR/$base"

  awk -F'\t' -v OFS='\t' '
    NR==1{
      for(i=1;i<=NF;i++){
        if(i==1 || tolower($i) !~ /putative/) keep[i]=1
      }
    }
    {
      first=1
      for(i=1;i<=NF;i++){
        if(keep[i]){
          if(first){ printf "%s",$i; first=0 }
          else     { printf "%s%s",OFS,$i }
        }
      }
      printf "\n"
    }
  ' "$f" > "$out"

  echo "Wrote: $out"
done

# quick check: confirm no 'putative' left in headers
grep -Hn "putative" "$OUTDIR"/*.tsv || echo "OK: no putative columns remain"
