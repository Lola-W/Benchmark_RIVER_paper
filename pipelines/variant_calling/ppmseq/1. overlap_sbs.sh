#!/usr/bin/env bash
# =============================================================================
# overlap_sbs.sh
#
# Purpose:
#   Keep only the outMap variants that are also present in sbsMap
#   (intersection of the two call sets).
#
# Usage:
#   bash overlap_sbs.sh <outMap.vcf.gz> <sbsMap.vcf.gz> <output.vcf.gz> [threads]
#
# Inputs:
#   outMap.vcf.gz   Bgzipped, indexed VCF of the primary call set
#   sbsMap.vcf.gz   Bgzipped, indexed VCF of the sbsMap call set
#
# Output:
#   output.vcf.gz   Records from outMap that are shared with sbsMap
#                   (bgzipped and tabix-indexed)
#
# Pipeline context:
#   Snakemake rule: overlap_sbs
#   The output is used downstream to train trinucleotide-context-specific
#   thresholds (see hom_and_single.py).
# =============================================================================

set -euo pipefail #Stop immediately on error

# Usage: bash overlap_sbs.sh <outMap.vcf.gz> <sbsMap.vcf.gz> <output.vcf.gz>

OUTMAP=$1
SBSMAP=$2
OUTPUT=$3
THREADS=${4:-10}  # Optional 4th argument; defaults to 10 threads

# -n=2 : keep records present in exactly both input files
# -w1  : write the records from the first file (outMap), so its original
#        INFO/FORMAT fields are preserved in the output
bcftools isec -n=2 -w1 -Oz --threads "$THREADS" \
  -o "$OUTPUT" \
  "$OUTMAP" "$SBSMAP"

# Index the output so downstream tools can use random access
tabix -f -p vcf "$OUTPUT" #indexing
