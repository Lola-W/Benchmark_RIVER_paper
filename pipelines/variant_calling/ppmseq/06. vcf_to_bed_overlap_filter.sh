#!/usr/bin/env bash
# =============================================================================
# vcf_to_bed_overlap_filter.sh
#
# Purpose:
#   Filter a VCF by overlap with a BED file. By default, keep variants that
#   overlap the BED regions; with -v, keep variants that do NOT overlap.
#
# Usage:
#   bash vcf_to_bed_overlap_filter.sh -i <input.vcf.gz> -b <regions.bed> \
#        -o <output.vcf.gz> [-v]
#
# Options:
#   -i  Input VCF (bgzipped)
#   -b  BED file of regions
#   -o  Output VCF (bgzipped; tabix-indexed automatically)
#   -v  Invert: keep NON-overlapping variants (default: keep overlapping)
#
# Notes:
#   - Overlap is evaluated at the variant's POS only (1-bp interval), not
#     across the full REF allele span of indels.
#   - Chromosome names must match between the VCF and the BED file.
#   - Requires bgzip, tabix, bedtools and awk.
#
# Pipeline context:
#   Snakemake rules step03_hmer (BED = homopolymers >= 7 bp) and
#   step04_str (BED = simple tandem repeats), both run with -v to
#   remove variants that fall in these regions.
# =============================================================================

set -euo pipefail

KEEP_NON_OVERLAP=0
IN="" BED="" OUT=""

while getopts ":i:b:o:v" opt; do
  case "$opt" in
    i) IN="$OPTARG" ;;
    b) BED="$OPTARG" ;;
    o) OUT="$OPTARG" ;;
    v) KEEP_NON_OVERLAP=1 ;;
  esac
done

# Work in a temporary directory that is removed on exit
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

HDR="$TMPDIR/header.vcf"
BODY="$TMPDIR/body.vcf"
BODYBED="$TMPDIR/body.bed"
PASSBED="$TMPDIR/passed.bed"

# Split the VCF into header lines and variant records in a single pass
bgzip -dc "$IN" | awk '
  /^#/ { print > "'"$HDR"'"; next }
         { print > "'"$BODY"'" }
'

# Convert each variant to a 1-bp BED interval (0-based start, 1-based end).
# The 4th column stores CHROM:POS:REF:ALT as a key to map results back to the VCF.
awk -F'\t' -v OFS='\t' '{
  print $1, $2-1, $2, $1":"$2":"$4":"$5
}' "$BODY" > "$BODYBED"

# Intersect variants with the regions
#   -v : report variants with NO overlap
#   -u : report each variant once if it overlaps at least one region
if [[ "$KEEP_NON_OVERLAP" -eq 1 ]]; then
  bedtools intersect -a "$BODYBED" -b "$BED" -v > "$PASSBED"
else
  bedtools intersect -a "$BODYBED" -b "$BED" -u > "$PASSBED"
fi

# Recover the original VCF records (all columns) for the retained keys
awk -F'\t' '
  NR==FNR { keep[$4]=1; next }
  { if ($1":"$2":"$4":"$5 in keep) print }
' "$PASSBED" "$BODY" > "$TMPDIR/filtered.body.vcf"

# Write the output: original header + retained records, then index
cat "$HDR" "$TMPDIR/filtered.body.vcf" | bgzip -c > "$OUT"
tabix -f -p vcf "$OUT"

# Report retained / total variant counts to the pipeline log
n_in=$(wc -l < "$BODY" | tr -d ' ')
n_out=$(wc -l < "$TMPDIR/filtered.body.vcf" | tr -d ' ')
[[ "$KEEP_NON_OVERLAP" -eq 1 ]] \
  && echo "[OK] Kept NON-overlapping variants: $n_out / $n_in" \
  || echo "[OK] Kept overlapping variants: $n_out / $n_in"
