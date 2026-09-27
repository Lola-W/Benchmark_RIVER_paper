#!/usr/bin/env bash
#SBATCH --job-name=scan2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=scan2-%j.log
set -euo pipefail
# Submit from pta/ with the archived SCAN2 environment activated.
mkdir -p output logs
for bam in input/*.sc.bam; do
  [[ -e "$bam" ]] || { echo 'No input/*.sc.bam files found' >&2; exit 1; }
  sample=$(basename "$bam" .sc.bam)
  mkdir -p "logs/$sample"
done
cd workflow
# Keep 24 threads per neuron as in the source workflow; process neurons sequentially.
snakemake --cores 24 --resources mem_mb=64000 --rerun-incomplete
