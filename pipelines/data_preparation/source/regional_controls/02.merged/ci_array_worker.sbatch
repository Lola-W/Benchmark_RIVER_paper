#!/usr/bin/env bash
#SBATCH --job-name=ci_worker
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00
#SBATCH --array=0-99
#SBATCH --export=ALL
#SBATCH --output=./work/regional_controls/02.merged/baseQ13_hg38/tmp/_ARRAY_LOGS/run_ci_array.%A_%a.log
#SBATCH --error=./work/regional_controls/02.merged/baseQ13_hg38/tmp/_ARRAY_LOGS/run_ci_array.%A_%a.err

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

# -----------------------------
# REQUIRED ENV (must be exported by submit script)
# -----------------------------
: "${CI_FASTA:?missing CI_FASTA}"
: "${CI_BAM:?missing CI_BAM}"
: "${CI_CI_PY:?missing CI_CI_PY}"
: "${CI_CHUNK_DIR:?missing CI_CHUNK_DIR}"
: "${CI_WORKDIR:?missing CI_WORKDIR}"
: "${CI_OUTDIR:?missing CI_OUTDIR}"

mkdir -p "$CI_WORKDIR" "$CI_OUTDIR"

IDX=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
VCF="${CI_CHUNK_DIR}/part_${IDX}.vcf"
MPILEUP_DIR="${CI_WORKDIR}/mpileup_${IDX}"
OUT_TSV="${CI_OUTDIR}/part_${IDX}.tsv"

mkdir -p "$MPILEUP_DIR"

python "$CI_CI_PY" "$CI_FASTA" "$CI_BAM" "$VCF" "$MPILEUP_DIR" "$OUT_TSV"

