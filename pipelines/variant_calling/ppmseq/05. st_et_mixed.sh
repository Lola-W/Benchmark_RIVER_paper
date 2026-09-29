#!/usr/bin/env bash
# =============================================================================
# st_et_mixed.sh
#
# Purpose:
#   Keep only variants whose strand-type (INFO/st) and error-type (INFO/et)
#   are both "MIXED".
#
# Usage:
#   bash st_et_mixed.sh <input.vcf.gz> <output.vcf.gz> [threads]
#
# Input:
#   input.vcf.gz    Bgzipped, indexed VCF with INFO/st and INFO/et tags
#
# Output:
#   output.vcf.gz   Variants with st == MIXED and et == MIXED
#                   (bgzipped and tabix-indexed)
#
# Pipeline context:
#   Snakemake rule: step02_st_et
# =============================================================================

set -euo pipefail

IN=$1
OUT=$2
THREADS=${3:-20}  # Optional 3rd argument; defaults to 20 threads

# Keep a record only if BOTH INFO/st and INFO/et equal "MIXED"
bcftools view \
  -i 'INFO/st=="MIXED" && INFO/et=="MIXED"' \
  -Oz --threads "$THREADS" \
  -o "$OUT" \
  "$IN"

# Index the output for downstream tools
tabix -f -p vcf "$OUT"

# Report the number of retained variants in the pipeline log
echo "[OK] wrote: $OUT"
bcftools index -n "$OUT" | awk '{print "[INFO] n_variants=" $1}'
