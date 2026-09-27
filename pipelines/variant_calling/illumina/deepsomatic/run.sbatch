#!/usr/bin/env bash
#SBATCH --job-name=deepsomatic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=deepsomatic-%j.log
set -euo pipefail
source "${1:?Pass config.sh as the first argument}"
mkdir -p "$OUTPUT_DIR/logs" "$OUTPUT_DIR/intermediate_results"
singularity run -B "$BIND_PATHS" "$IMAGE_PATH" run_deepsomatic \
  --model_type=WGS --ref="$REF_FASTA" \
  --reads_normal="$NORMAL_BAM_PATH" --reads_tumor="$TUMOR_BAM_PATH" \
  --output_vcf="$OUTPUT_DIR/7614_300x_cpu_ver.vcf.gz" \
  --output_gvcf="$OUTPUT_DIR/7614_300x_cpu_ver.g.vcf.gz" \
  --sample_name_tumor="$(basename "$TUMOR_BAM_PATH" .bam)" \
  --sample_name_normal="$(basename "$NORMAL_BAM_PATH" .bam)" \
  --num_shards="${SLURM_CPUS_PER_TASK:-$(nproc)}" \
  --logging_dir="$OUTPUT_DIR/logs" \
  --intermediate_results_dir="$OUTPUT_DIR/intermediate_results" \
  --use_default_pon_filtering=false --dry_run=false
