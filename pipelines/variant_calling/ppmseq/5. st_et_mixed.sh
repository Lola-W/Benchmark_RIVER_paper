#!/usr/bin/env bash
set -euo pipefail

# Usage: bash st_et_mixed.sh <input.vcf.gz> <output.vcf.gz>

IN=$1
OUT=$2
THREADS=${3:-20}


# Keep only variants where both strand type(st) and error type(et) are MIXED
bcftools view \
  -i 'INFO/st=="MIXED" && INFO/et=="MIXED"' \
  -Oz --threads "$THREADS" \
  -o "$OUT" \
  "$IN"

tabix -f -p vcf "$OUT"

echo "[OK] wrote: $OUT"
bcftools index -n "$OUT" | awk '{print "[INFO] n_variants=" $1}'
