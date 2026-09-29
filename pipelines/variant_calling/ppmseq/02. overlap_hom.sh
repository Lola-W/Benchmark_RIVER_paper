#!/usr/bin/env bash
# =============================================================================
# overlap_hom.sh
#
# Purpose:
#   Keep only the outMap variants that are also present in homMap
#   (intersection of the two call sets).
#
# Usage:
#   bash overlap_hom.sh <outMap.vcf.gz> <homMap.vcf.gz> <output.vcf.gz> [threads]
#
# Inputs:
#   outMap.vcf.gz   Bgzipped, indexed VCF of the primary call set
#   homMap.vcf.gz   Bgzipped, indexed VCF of the homMap call set
#
# Output:
#   output.vcf.gz   Records from outMap that are shared with homMap
#                   (bgzipped and tabix-indexed)
#
# Pipeline context:
#   Snakemake rule: overlap_hom
#   The output is used downstream to train trinucleotide-context-specific
#   thresholds (see hom_and_single.py).
# =============================================================================

set -euo pipefail


OUTMAP=$1
HOMMAP=$2
OUTPUT=$3
THREADS=${4:-10}  # Optional 4th argument; defaults to 10 threads

# -n=2 : keep records present in exactly both input files
# -w1  : write the records from the first file (outMap), so its original
#        INFO/FORMAT fields are preserved in the output
bcftools isec -n=2 -w1 -Oz --threads "$THREADS" \
  -o "$OUTPUT" \
  "$OUTMAP" "$HOMMAP"

tabix -f -p vcf "$OUTPUT"
