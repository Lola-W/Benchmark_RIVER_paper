#!/usr/bin/env bash
#SBATCH --job-name=dupcaller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --output=dupcaller-%j.log
set -euo pipefail
source "${1:?Pass nano.config.sh or ud.config.sh}"
mkdir -p "$(dirname "$OUTPREFIX")"
DupCaller.py call -b "$TUMOR" -f "$REF" -o "$OUTPREFIX" \
  -p "${SLURM_CPUS_PER_TASK:-24}" -n "$NORMAL" -g "$GERMLINE" -m "$NOISE"
