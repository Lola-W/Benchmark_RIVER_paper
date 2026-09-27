#!/usr/bin/env bash
#SBATCH --job-name=concat_ibp
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --export=ALL

set -euo pipefail

mkdir -p all_ibp

tmp="all_ibp/${METHOD}_ibp.tsv.tmp"
out="all_ibp/${METHOD}_ibp.tsv"

cat "${SRCDIR}"/chr{1..22}.tsv \
    "${SRCDIR}"/chrX.tsv \
    > "$tmp"

mv "$tmp" "$out"

echo "[DONE] ${METHOD}"
