#!/usr/bin/env bash
#SBATCH --job-name=duplex_clean_chr
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00
#SBATCH --array=0-22
#SBATCH --export=ALL
#SBATCH --output=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/duplex_clean.%A_%a.log
#SBATCH --error=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/duplex_clean.%A_%a.err

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

# -----------------------------
# REQUIRED ENV (export these when submitting)
# -----------------------------
: "${IN_BAM:?missing IN_BAM}"
: "${OUT_BASEDIR:?missing OUT_BASEDIR}"
: "${PY_SCRIPT:?missing PY_SCRIPT}"

mkdir -p "${OUT_BASEDIR}/tmp/_CHR_LOGS"
mkdir -p "${OUT_BASEDIR}/tmp/by_chr"

BAM_BASENAME="$(basename "${IN_BAM}")"
BAM_STEM="${BAM_BASENAME%.bam}"

CHRS=(chr{1..22} chrX)
CHR="${CHRS[$SLURM_ARRAY_TASK_ID]}"

WORKDIR="${OUT_BASEDIR}/tmp/by_chr/${CHR}"
mkdir -p "${WORKDIR}"

# Unsorted from python, then sort+index here
OUT_BAM_UNSORTED="${WORKDIR}/duplex_clean.${BAM_STEM}.${CHR}.unsorted.bam"
OUT_BAM="${OUT_BASEDIR}/duplex_clean.${BAM_STEM}.${CHR}.bam"

# -----------------------------
# 0) Sanity check only on task 0 (chr1)
# -----------------------------
if [[ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]]; then
  SANITY_BAM="${WORKDIR}/sanity.duplex_clean.${BAM_STEM}.chr1_1999800_2000500.unsorted.bam"
  python3 "${PY_SCRIPT}" \
    --bam "${IN_BAM}" \
    --region "chr1:1999800-2000500" \
    --out-bam "${SANITY_BAM}" \
    --threads 4 \
    --progress-every 1000000

  # sort+index sanity
  samtools sort -@ 4 -m 4G -o "${SANITY_BAM%.unsorted.bam}.bam" "${SANITY_BAM}"
  samtools index -@ 4 "${SANITY_BAM%.unsorted.bam}.bam"
  rm -f "${SANITY_BAM}"
fi

# -----------------------------
# 1) Per-chrom duplex clean
# -----------------------------
python3 "${PY_SCRIPT}" \
  --bam "${IN_BAM}" \
  --chrom "${CHR}" \
  --out-bam "${OUT_BAM_UNSORTED}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --progress-every 5000000

# -----------------------------
# 2) Sort+index (so final concat is safe/fast)
# -----------------------------
TMP_SORT="${WORKDIR}/duplex_clean.${BAM_STEM}.${CHR}.sorted.bam"
samtools sort -@ "${SLURM_CPUS_PER_TASK}" -m 6G -o "${TMP_SORT}" "${OUT_BAM_UNSORTED}"
mv "${TMP_SORT}" "${OUT_BAM}"
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${OUT_BAM}"

rm -f "${OUT_BAM_UNSORTED}"

echo "[DONE] ${CHR} -> ${OUT_BAM} (+.bai)"
