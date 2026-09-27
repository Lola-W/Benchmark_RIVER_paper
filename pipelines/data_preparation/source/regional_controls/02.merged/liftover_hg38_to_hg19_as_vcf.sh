#!/usr/bin/env bash
set -euo pipefail

# usage:
#   bash liftover_hg38_to_hg19_as_vcf.sh big_with_downsampling_ppmSeq_extended.annot.tsv
#
# output:
#   <input>.with_hg19.tsv
# columns appended:
#   hg19_chr hg19_pos hg19_ref hg19_alt hg19_unmapped

IN="${1:-}"
CHAIN="${CHAIN:-hg38ToHg19.over.chain.gz}"
REF37="${REF37:-./work/regional_controls/human_g1k_v37_decoy.fasta}"
KEEP_TMP="${KEEP_TMP:-0}"

if [[ -z "$IN" ]]; then
  echo "Usage: $0 <input.tsv>" >&2
  exit 1
fi
if [[ ! -f "$IN" ]]; then
  echo "[ERROR] input not found: $IN" >&2
  exit 1
fi
if [[ ! -f "$CHAIN" ]]; then
  echo "[ERROR] chain not found: $CHAIN" >&2
  exit 1
fi
if [[ ! -f "$REF37" ]]; then
  echo "[ERROR] hg19 ref not found: $REF37" >&2
  exit 1
fi

OUT="${IN%.tsv}.with_hg19.tsv"

tmpdir="$(mktemp -d)"
if [[ "$KEEP_TMP" == "1" ]]; then
  echo "[INFO] KEEP_TMP=1 tmpdir=$tmpdir" >&2
else
  trap 'rm -rf "$tmpdir"' EXIT
fi

vcf_in="$tmpdir/in.vcf"
vcf_out="$tmpdir/out.vcf"
map="$tmpdir/id2hg19.tsv"
crosslog="$tmpdir/crossmap.log"

# -------------------------
# 1) TSV -> fake hg38 VCF (only VCF-valid alleles)
#    - Keep row id (RID) = body row number (NR)
# -------------------------
{
  echo "##fileformat=VCFv4.2"
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  tail -n +2 "$IN" \
  | awk -v OFS="\t" '
    BEGIN{ IGNORECASE=1 }
    {
      id=NR
      chr=$1; pos=$2; ref=$3; alt=$4

      # POS must be positive integer
      if (pos !~ /^[0-9]+$/ || pos < 1) next

      # VCF alleles must be non-empty and not "."
      if (ref=="." || alt=="." || ref=="" || alt=="") next

      # Restrict to A/C/G/T/N only (avoid weird symbols breaking CrossMap)
      if (ref !~ /^[ACGTN]+$/ || alt !~ /^[ACGTN]+$/) next

      print chr, pos, id, ref, alt, ".", "PASS", "RID="id
    }
  '
} > "$vcf_in"

# -------------------------
# 2) CrossMap vcf hg38->hg19 (do NOT swallow errors)
# -------------------------
CrossMap vcf "$CHAIN" "$vcf_in" "$REF37" "$vcf_out" 2> "$crosslog" || {
  echo "[ERROR] CrossMap failed. Last 200 lines of log:" >&2
  tail -n 200 "$crosslog" >&2
  exit 1
}

# output might be out.vcf or out.vcf.gz (if CrossMap decided to compress)
VCF_MAPPED=""
if [[ -s "$vcf_out" ]]; then
  VCF_MAPPED="$vcf_out"
elif [[ -s "${vcf_out}.gz" ]]; then
  VCF_MAPPED="${vcf_out}.gz"
else
  echo "[ERROR] CrossMap produced no output VCF: $vcf_out(.gz)" >&2
  echo "[INFO] tmpdir contents:" >&2
  ls -lh "$tmpdir" >&2 || true
  exit 1
fi

# -------------------------
# 3) mapped VCF -> id map: id -> hg19 chr,pos,ref,alt
#    - hg19 chr: 1-22,X,Y,MT (strip chr, M->MT)
# -------------------------
if [[ "$VCF_MAPPED" == *.gz ]]; then
  zcat "$VCF_MAPPED" \
  | awk -v OFS="\t" '
      BEGIN{ IGNORECASE=1 }
      !/^#/ {
        id=$3
        chr=$1
        sub(/^chr/,"",chr)
        if (chr=="M") chr="MT"
        print id, chr, $2, $4, $5
      }
    ' > "$map"
else
  awk -v OFS="\t" '
      BEGIN{ IGNORECASE=1 }
      !/^#/ {
        id=$3
        chr=$1
        sub(/^chr/,"",chr)
        if (chr=="M") chr="MT"
        print id, chr, $2, $4, $5
      }
    ' "$VCF_MAPPED" > "$map"
fi

# -------------------------
# 4) merge back to TSV (preserve order; do not drop)
# -------------------------
awk -v OFS="\t" -v MAP="$map" '
  BEGIN {
    while ((getline < MAP) > 0)
      m[$1]=$2"\t"$3"\t"$4"\t"$5
    close(MAP)
  }
  NR==1 {
    print $0, "hg19_chr","hg19_pos","hg19_ref","hg19_alt","hg19_unmapped"
    next
  }
  {
    id=NR-1
    if (id in m) {
      split(m[id],a,"\t")
      print $0,a[1],a[2],a[3],a[4],0
    } else {
      print $0,".",".",".",".",1
    }
  }
' "$IN" > "$OUT"

echo "[OK] $OUT"

