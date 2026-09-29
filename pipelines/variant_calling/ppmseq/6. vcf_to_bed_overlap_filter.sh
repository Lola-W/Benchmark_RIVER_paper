#!/usr/bin/env bash
set -euo pipefail

# Usage: bash vcf_to_bed_overlap_filter.sh -i <input.vcf.gz> -b <regions.bed> -o <output.vcf.gz> [-v]
# -v: keep NON-overlapping variants (default: keep overlapping)

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

# Split header and body in a single pass
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

HDR="$TMPDIR/header.vcf"
BODY="$TMPDIR/body.vcf"
BODYBED="$TMPDIR/body.bed"
PASSBED="$TMPDIR/passed.bed"

bgzip -dc "$IN" | awk '
  /^#/ { print > "'"$HDR"'"; next }
         { print > "'"$BODY"'" }
'

# Convert VCF body to 1bp BED with key (CHROM:POS:REF:ALT)
awk -F'\t' -v OFS='\t' '{
  print $1, $2-1, $2, $1":"$2":"$4":"$5
}' "$BODY" > "$BODYBED"

# Intersect with BED
if [[ "$KEEP_NON_OVERLAP" -eq 1 ]]; then
  bedtools intersect -a "$BODYBED" -b "$BED" -v > "$PASSBED"
else
  bedtools intersect -a "$BODYBED" -b "$BED" -u > "$PASSBED"
fi

# Map back to original VCF records
awk -F'\t' '
  NR==FNR { keep[$4]=1; next }
  { if ($1":"$2":"$4":"$5 in keep) print }
' "$PASSBED" "$BODY" > "$TMPDIR/filtered.body.vcf"

# Write output
cat "$HDR" "$TMPDIR/filtered.body.vcf" | bgzip -c > "$OUT"
tabix -f -p vcf "$OUT"

# Report
n_in=$(wc -l < "$BODY" | tr -d ' ')
n_out=$(wc -l < "$TMPDIR/filtered.body.vcf" | tr -d ' ')
[[ "$KEEP_NON_OVERLAP" -eq 1 ]] \
  && echo "[OK] Kept NON-overlapping variants: $n_out / $n_in" \
  || echo "[OK] Kept overlapping variants: $n_out / $n_in"
