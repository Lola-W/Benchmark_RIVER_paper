#!/usr/bin/env bash
# =============================================================================
# make_final_vcfs.sh
#
# Purpose:
#   Produce the three final ppmSeq call sets from the multiread, singleton
#   and HC-singleton VCFs of trinucDenoised_integration_refining.py:
#     - final_multiread            multiread variants, non-multiallelic in raw data
#     - final_putative_multiread   singleton after filtering, but seen >= 2 times
#                                  in the raw outMap (non-multiallelic)
#     - final_singleton_HC         singletons seen exactly once at every level
#
# Usage:
#   bash make_final_vcfs.sh <sample_name> <multiread.vcf.gz> <singleton.vcf.gz> \
#        <singleton_hc.vcf.gz> <output_dir>
#
# Inputs:
#   multiread.vcf.gz     DUP_COUNT_FILTERED >= 2   (step07-1)
#   singleton.vcf.gz     DUP_COUNT_FILTERED == 1   (step07-2)
#   singleton_hc.vcf.gz  candidate HC singletons   (step08-2)
#   All carry the INFO/DUP_COUNT_* annotations added in the previous step.
#
# Outputs (in <output_dir>, bgzipped and tabix-indexed):
#   final_multiread_<sample>.vcf.gz
#   final_putative_multiread_<sample>.vcf.gz
#   final_singleton_HC_<sample>.vcf.gz
#
# Notes:
#   - "Non-multiallelic" means DUP_COUNT_RAW_MULTI_ALLELE == DUP_COUNT_RAW,
#     i.e. every raw outMap record at that position carries the same
#     REF/ALT as the variant.
#   - Per-filter before/after counts are printed to the pipeline log.
#
# Pipeline context:
#   Snakemake rule: final (last step).
# =============================================================================

set -euo pipefail


SAMPLENAME=$1
MULTIREAD_IN=$2
SINGLETON_IN=$3
SINGLETON_HC_IN=$4
OUTDIR=$5

FINAL_MULTIREAD="${OUTDIR}/final_multiread_${SAMPLENAME}.vcf.gz"
FINAL_PUTATIVE_MULTIREAD="${OUTDIR}/final_putative_multiread_${SAMPLENAME}.vcf.gz"
FINAL_SINGLETON_HC="${OUTDIR}/final_singleton_HC_${SAMPLENAME}.vcf.gz"

# Helpers: count records, index a VCF, print a before/after summary line
count() { bcftools view -H "$1" | wc -l; }
idx()   { tabix -f -p vcf "$1"; }
report() {
    local before=$1 after=$2 label=$3
    echo "[$label] Before=$before / After=$after / Removed=$((before - after))"
}

echo "=== $SAMPLENAME ==="

# 1. Multiread: keep non-multiallelic only (DUP_COUNT_RAW_MULTI_ALLELE == DUP_COUNT_RAW)
before=$(count "$MULTIREAD_IN")
bcftools view -i 'INFO/DUP_COUNT_RAW_MULTI_ALLELE=INFO/DUP_COUNT_RAW' "$MULTIREAD_IN" -Oz -o "$FINAL_MULTIREAD"
idx "$FINAL_MULTIREAD"
report "$before" "$(count "$FINAL_MULTIREAD")" "final_multiread"

# 2. Putative multiread: singletons after filtering that appear >= 2 times in the raw outMap
before=$(count "$SINGLETON_IN")
bcftools view -i 'INFO/DUP_COUNT_RAW>=2 && INFO/DUP_COUNT_RAW_MULTI_ALLELE=INFO/DUP_COUNT_RAW' "$SINGLETON_IN" -Oz -o "$FINAL_PUTATIVE_MULTIREAD"
idx "$FINAL_PUTATIVE_MULTIREAD"
report "$before" "$(count "$FINAL_PUTATIVE_MULTIREAD")" "final_putative_multiread"

# 3. Singleton HC: QC check then keep only exact 1,1,1,1 ((DUP_COUNT_ 1. FILTERED, 2. RAW_MULTI_ALLELE, 3. RAW, 4. SNVQ40)
echo "[singleton_HC QC]"
bcftools query -f '%INFO/DUP_COUNT_FILTERED\t%INFO/DUP_COUNT_RAW_MULTI_ALLELE\t%INFO/DUP_COUNT_RAW\t%INFO/DUP_COUNT_SNVQ40\n' "$SINGLETON_HC_IN" \
| awk -F'\t' '
    { dcf=$1+0; raw_multi=$2+0; raw=$3+0; snvq40=$4+0
      if (dcf==1 && raw_multi==1 && raw==1 && snvq40==1) all1111++
      else not1111++
      if (raw_multi != raw) multiallelic++ }
    END { print "all 1,1,1,1:", all1111+0
          print "not all 1,1,1,1:", not1111+0
          print "raw-level multiallelic:", multiallelic+0 }'

before=$(count "$SINGLETON_HC_IN")
bcftools view -i 'INFO/DUP_COUNT_FILTERED=1 && INFO/DUP_COUNT_RAW_MULTI_ALLELE=1 && INFO/DUP_COUNT_RAW=1 && INFO/DUP_COUNT_SNVQ40=1' "$SINGLETON_HC_IN" -Oz -o "$FINAL_SINGLETON_HC"
idx "$FINAL_SINGLETON_HC"
report "$before" "$(count "$FINAL_SINGLETON_HC")" "final_singleton_HC"

echo ""
echo "[final files]"
ls -lh "$FINAL_MULTIREAD" "$FINAL_MULTIREAD.tbi" \
       "$FINAL_PUTATIVE_MULTIREAD" "$FINAL_PUTATIVE_MULTIREAD.tbi" \
       "$FINAL_SINGLETON_HC" "$FINAL_SINGLETON_HC.tbi"
