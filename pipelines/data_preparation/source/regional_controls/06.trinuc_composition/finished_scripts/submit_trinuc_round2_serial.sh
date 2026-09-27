#!/usr/bin/env bash
set -euo pipefail

BASE="./work/regional_controls"
FA="${BASE}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BLACK="${BASE}/06.trinuc_composition/all_blacklist_hg38.bed"
OUTROOT="${BASE}/06.trinuc_composition/baseQ13_hg38_v4"

WORKER="trinuc_round2_serial.sbatch"

declare -A BAM

CLEAN_DIR="${BASE}/05.cleaned_bam"

# NanoSeq
BAM[nanoseq_clean_7614]="${CLEAN_DIR}/duplex_clean.NanoSeq_7614.chr1-22_X.bam"
BAM[nanoseq_clean_6566]="${CLEAN_DIR}/duplex_clean.NanoSeq_6566.chr1-22_X.bam"

# UDSeq
BAM[udseq_clean_7614]="${CLEAN_DIR}/duplex_clean.UDSeq_7614.chr1-22_X.bam"
BAM[udseq_clean_30x]="${CLEAN_DIR}/duplex_clean.30x_UDSeq_7614.chr1-22_X.bam"
BAM[udseq_clean_6566]="${CLEAN_DIR}/duplex_clean.UDSeq_6566.chr1-22_X.bam"

# ppmSeq
BAM[ppmseq_7614]="${CLEAN_DIR}/ppmSeq_7614.mixed_mixed.bam"
BAM[ppmseq_30x]="${CLEAN_DIR}/ppmSeq_30x_7614.mixed_mixed.bam"
BAM[ppmseq_100x]="${CLEAN_DIR}/ppmSeq_100x_7614.mixed_mixed.bam"
BAM[ppmseq_200x]="${CLEAN_DIR}/ppmSeq_200x_7614.mixed_mixed.bam"
BAM[ppmseq_6566]="${CLEAN_DIR}/ppmSeq_6566.mixed_mixed.bam"

# HiDEF
BAM[hidef_clean_7614]="${CLEAN_DIR}/hidef_duplex.7614.7614.7614.ccs.filtered.aligned.sorted.chr1-22_X.bam"

# Illumina
BAM[illumina_full]="./work/project/data/illumina_hg38_bam/7614_illumina_brain_somatic.bam"
BAM[illumina_30x]="./work/project/data/illumina_hg38_bam/illumina_30x.bam"
BAM[illumina_100x]="./work/project/data/illumina_hg38_bam/illumina_100x.bam"
BAM[illumina_200x]="./work/project/data/illumina_hg38_bam/illumina_200x.bam"
BAM[illumina_6566_blood]="${BASE}/00.bam_hg38/6566_blood_hg38.bam"
BAM[illumina_6566_sperm]="${BASE}/00.bam_hg38/6566_sperm_hg38.bam"

# PTA
PTA_DIR="./inputs/pta_hg38"
for i in {1..16}; do
  key="PTA_$(printf "%02d" $i)"
  BAM["$key"]="${PTA_DIR}/7614_single_neuron_${i}.bam"
done

for key in "${!BAM[@]}"; do
  echo "[submit round2 serial] $key"
  sbatch --job-name="trinuc_r2_${key}" \
    --export=FA="$FA",BLACK="$BLACK",BAM="${BAM[$key]}",KEY="$key",OUTROOT="$OUTROOT" \
    "$WORKER"
done
