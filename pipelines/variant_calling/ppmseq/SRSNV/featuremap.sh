#!/bin/bash
# featuremap.sh
# Create a FeatureMap for one CRAM (or one genomic shard of it) with
# GATK FlowFeatureMapper, then apply VariantFiltration.
#
# Usage:
#   bash featuremap.sh <sample> <shard> <input.cram> <out_base> \
#        [intervals] [padding] [mem_gb] [threads]
#
# Environment variables:
#   REF        hg38 reference FASTA
#   GATK_JAR   GATK jar (tested with GATK 4.6.2.0)
#   JAVA       java executable (default: java)
set -euo pipefail

: "${REF:?Set REF to the reference FASTA}"
: "${GATK_JAR:?Set GATK_JAR to the GATK jar}"
JAVA="${JAVA:-java}"
MEM=9

trim() { sed -e 's/^[[:space:]]\+//' -e 's/[[:space:]]\+$//'; }

cramBase="$(printf '%s' "${1:-}" | trim)"
shard="$(printf '%s' "${2:-}" | trim)"
input_cram="$(printf '%s' "${3:-}" | trim)"
OUT_BASE="$(printf '%s' "${4:-}" | trim)"
INTERVALS="$(printf '%s' "${5:-}" | trim)"
PAD="$(printf '%s' "${6:-0}" | trim)"
MEM_GB="$(printf '%s' "${7:-$MEM}" | trim)"
THREADS="$(printf '%s' "${8:-16}" | trim)"
export OMP_NUM_THREADS="${THREADS}"

SAMPLE_DIR="${OUT_BASE}/18.srsnv/${cramBase}"
FEATURE_DIR="${SAMPLE_DIR}/00.featureMap"
mkdir -p "${FEATURE_DIR}"

tmp_feature_vcf="${FEATURE_DIR}/${cramBase}_${shard}.fmm.tmp.vcf.gz"
feature_vcf="${FEATURE_DIR}/${cramBase}_${shard}.feature.vcf.gz"

# Restrict to the shard's intervals if given
FMM_L=""
VF_L=""
if [[ -n "${INTERVALS}" ]]; then
  FMM_L="-L ${INTERVALS} -ip ${PAD}"
  VF_L="-L ${INTERVALS}"
fi

echo "[$(date +'%F %T')] FlowFeatureMapper start sample=${cramBase} shard=${shard} mem=${MEM_GB}G threads=${THREADS}"
"${JAVA}" -Xmx${MEM_GB}g -jar "${GATK_JAR}" FlowFeatureMapper \
  -I "${input_cram}" -O "${tmp_feature_vcf}" -R "${REF}" ${FMM_L} \
  --snv-identical-bases 5 --snv-identical-bases-after 5 --min-score 0 --limit-score 10 \
  --read-filter MappingQualityReadFilter --minimum-mapping-quality 60 \
  --flow-use-t0-tag --flow-fill-empty-bins-value 0.0001 --surrounding-median-quality-size 20 \
  --copy-attr tm --copy-attr a3 --copy-attr rq --copy-attr st --copy-attr et \
  --copy-attr as --copy-attr ts --copy-attr ae --copy-attr te --copy-attr s3 --copy-attr s2

echo "[$(date +'%F %T')] VariantFiltration start sample=${cramBase} shard=${shard}"
"${JAVA}" -Xmx${MEM_GB}g -jar "${GATK_JAR}" VariantFiltration \
  -R "${REF}" -V "${tmp_feature_vcf}" -O "${feature_vcf}" ${VF_L} \
  --create-output-variant-index true

rm -f "${tmp_feature_vcf}" "${tmp_feature_vcf}.tbi" 2>/dev/null || true
echo "[$(date +'%F %T')] Done -> ${feature_vcf}"
