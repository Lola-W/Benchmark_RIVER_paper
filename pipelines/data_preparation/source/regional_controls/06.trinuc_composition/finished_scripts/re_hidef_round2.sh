#!/usr/bin/env bash
set -euo pipefail

BASE="./work/regional_controls"
FA="${BASE}/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BLACK="${BASE}/06.trinuc_composition/all_blacklist_hg38.bed"
OUTROOT="${BASE}/06.trinuc_composition/baseQ13_hg38_v4"

WORKER="${BASE}/06.trinuc_composition/trinuc_round2_serial.sbatch"

PTA_DIR="./inputs/pta_hg38"


key=hidef_rerun
bam="./work/regional_controls/05.cleaned_bam/7614.7614.7614.ccs.filtered.aligned.sorted.bam"

echo "[rerun round2] $key"

sbatch --job-name="trinuc_r2_${key}_rerun" \
	--export=FA="$FA",BLACK="$BLACK",BAM="$bam",KEY="$key",OUTROOT="$OUTROOT" "$WORKER"
