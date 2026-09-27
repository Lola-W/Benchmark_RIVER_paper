#!/usr/bin/env bash
#SBATCH --job-name=deepmosaic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=84G
#SBATCH --time=7-00:00:00
#SBATCH --output=deepmosaic-%j.log
set -euo pipefail
source "${1:?Pass config.sh as the first argument}"
INPUT_TSV="$(python3 "$PACKAGE_ROOT/prepare_paths.py" "$INPUT_TSV")"
mkdir -p "$OUTPUT_DIR"
singularity exec -B "$BIND_PATHS" "$IMAGE_PATH" deepmosaic-draw \
  -i "$INPUT_TSV" -o "$OUTPUT_DIR" -a "$ANNOVAR_DIR" -b hg38 -db refGene
singularity exec -B "$BIND_PATHS" "$IMAGE_PATH" deepmosaic-predict \
  -i "$OUTPUT_DIR/features.txt" -o "$OUTPUT_DIR/output.txt" \
  -m efficientnet-b4_epoch_6.pt -b 12 -gb hg38
