#!/usr/bin/env bash
set -euo pipefail

BASE="./work/regional_controls"
FA="${BASE}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BLACK="${BASE}/06.trinuc_composition/all_blacklist_hg38.bed"
OUTROOT="${BASE}/06.trinuc_composition/baseQ13_hg38_v4"

WORKER="${BASE}/06.trinuc_composition/trinuc_round2_serial.sbatch"

PTA_DIR="./inputs/pta_hg38"

for i in 12 14 15; do
  key="PTA_$(printf "%02d" "$i")"
  bam="${PTA_DIR}/7614_single_neuron_${i}.bam"

  # clean previous tmp/output
  rm -rf "${OUTROOT}/round2_tmp/${key}"
  rm -f  "${OUTROOT}/round2_allchr/${key}.allchr.tsv"

  echo "[rerun round2] $key"

  sbatch --job-name="trinuc_r2_${key}_rerun" \
    --export=FA="$FA",BLACK="$BLACK",BAM="$bam",KEY="$key",OUTROOT="$OUTROOT" \
    "$WORKER"
done
