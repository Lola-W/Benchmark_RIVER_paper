#!/usr/bin/env bash
#SBATCH --job-name=ref_tri32
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-22
#SBATCH --output=log/%x.%A_%a.out
#SBATCH --error=log/%x.%A_%a.err
#SBATCH --export=ALL

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

export LC_ALL=C

BASE="./work/regional_controls/06.trinuc_composition"

FASTA="./work/regional_controls/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BED="${BASE}/all_blacklist_hg38.bed"

OUTDIR="${BASE}/01.trinuc32/by_chrom_tmp"
mkdir -p "${OUTDIR}" log

CHR_LIST=(chr{1..22} chrX)
CHR="${CHR_LIST[${SLURM_ARRAY_TASK_ID}]}"

OUT_TMP="${OUTDIR}/ref_hg38_trinuc.${CHR}.tsv.tmp"
OUT="${OUTDIR}/ref_hg38_trinuc.${CHR}.tsv"

python "${BASE}/ref_trinuc32_onechr.py" \
  --fasta "${FASTA}" \
  --bed "${BED}" \
  --chrom "${CHR}" \
  --out "${OUT_TMP}"

mv "${OUT_TMP}" "${OUT}"
echo "[DONE] ${CHR} -> ${OUT}"
