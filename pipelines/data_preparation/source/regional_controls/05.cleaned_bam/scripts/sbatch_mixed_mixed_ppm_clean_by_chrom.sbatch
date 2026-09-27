#!/usr/bin/env bash
#SBATCH --job-name=ppm_mixmix_chr
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
#SBATCH --output=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/ppm_mixmix.%A_%a.log
#SBATCH --error=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/ppm_mixmix.%A_%a.err

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

: "${IN_CRAM:?missing IN_CRAM}"
: "${OUT_BASEDIR:?missing OUT_BASEDIR}"

mkdir -p "${OUT_BASEDIR}/tmp/_CHR_LOGS"
mkdir -p "${OUT_BASEDIR}/tmp/by_chr"

BASENAME="$(basename "${IN_CRAM}")"
STEM="${BASENAME%.cram}"
STEM="${STEM%.bam}"

CHRS=(chr{1..22} chrX)
CHR="${CHRS[$SLURM_ARRAY_TASK_ID]}"

WORKDIR="${OUT_BASEDIR}/tmp/by_chr/${STEM}/${CHR}"
mkdir -p "${WORKDIR}"

EXCL_FLAGS=3840
EXPR='[st]=="MIXED" && [et]=="MIXED"'

OUT_UNSORTED="${WORKDIR}/ppm_mixedmixed.${STEM}.${CHR}.unsorted.bam"
OUT_BAM="${OUT_BASEDIR}/ppm_mixedmixed.${STEM}.${CHR}.bam"

echo "Processing ${STEM} ${CHR}"

# 1. Extract MIXED-MIXED reads
samtools view -h \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -F "${EXCL_FLAGS}" \
  -e "${EXPR}" \
  "${IN_CRAM}" "${CHR}" \
| samtools view \
    -@ "${SLURM_CPUS_PER_TASK}" \
    -b -o "${OUT_UNSORTED}" -

# 2. Sort
TMP_SORT="${WORKDIR}/sorted.bam"
samtools sort \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -m 6G \
  -o "${TMP_SORT}" \
  "${OUT_UNSORTED}"

mv "${TMP_SORT}" "${OUT_BAM}"
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${OUT_BAM}"

rm -f "${OUT_UNSORTED}"

echo "[DONE] ${STEM} ${CHR}"
