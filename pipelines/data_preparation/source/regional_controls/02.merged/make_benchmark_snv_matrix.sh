#!/usr/bin/env bash
set -euo pipefail

BIG="big_with_downsampling_ppmSeq_extended.annot.tsv"
ILL="snv_illumina_detected_presence_by_depth.tsv"
OUT="benchmark_snv_matrix.tsv"

awk -v FS="\t" -v OFS="\t" '
function key(c,p,r,a){ return c "\t" p "\t" r "\t" a }

# 1) load Illumina presence table
NR==FNR {
  if (FNR==1) {
    for (i=1;i<=NF;i++) H[$i]=i
    next
  }
  k = key($H["chr"], $H["pos"], $H["ref"], $H["alt"])
  ILL[k] = $(H["Illumina"]) OFS $(H["Illumina_30x"]) OFS $(H["Illumina_100x"]) OFS $(H["Illumina_200x"])
  next
}

# 2) big header: find crop column index up to AmpliSeq
FNR==1 {
  end_col = 0
  for (i=1;i<=NF;i++) if ($i=="AmpliSeq") end_col=i
  if (end_col==0) end_col=NF

  # print cropped header
  for (i=1;i<=end_col;i++) {
    printf "%s%s", $i, (i==end_col?ORS:OFS)
  }
  next
}

# 3) big data: overwrite illumina cols (5-8) from ILL; if absent, set to 0
{
  k = key($1,$2,$3,$4)
  seen_big[k]=1

  if (k in ILL) {
    split(ILL[k], v, OFS)
    $5=v[1]; $6=v[2]; $7=v[3]; $8=v[4]
  } else {
    $5=0; $6=0; $7=0; $8=0
  }

  for (i=1;i<=end_col;i++) {
    printf "%s%s", $i, (i==end_col?ORS:OFS)
  }
}

# 4) add Illumina-only variants (new rows)
END {
  for (k in ILL) {
    if (!(k in seen_big)) {
      split(k, f, "\t")
      split(ILL[k], v, OFS)

      # columns 1-4
      printf "%s\t%s\t%s\t%s", f[1],f[2],f[3],f[4]
      # columns 5-8 (Illumina updated)
      printf "\t%s\t%s\t%s\t%s", v[1],v[2],v[3],v[4]
      # columns 9..end_col (all zeros)
      for (i=9;i<=end_col;i++) printf "\t0"
      printf "\n"
    }
  }
}
' "$ILL" "$BIG" > "$OUT"

