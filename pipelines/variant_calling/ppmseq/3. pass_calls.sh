#!/usr/bin/env bash
# =============================================================================
# pass_calls.sh
#
# Purpose:
#   Extract variants with FILTER == PASS from the outMap VCF.
#
# Usage:
#   bash pass_calls.sh <outMap.vcf.gz> <output.vcf.gz>
#
# Input:
#   outMap.vcf.gz   Bgzipped VCF of the primary call set
#
# Output:
#   output.vcf.gz   PASS-only variants (bgzipped and tabix-indexed)
#
# Pipeline context:
#   Snakemake rule: step00_pass (first filtering step)
# =============================================================================

set -euo pipefail


IN=$1 #input: outMap.vcf.gz
OUT=$2 #output: PASS-only variants

# -f PASS : keep only records whose FILTER field is PASS
bcftools view -f PASS "$IN" -Oz -o "$OUT"

# Index the output for downstream tools
tabix -f -p vcf "$OUT"
