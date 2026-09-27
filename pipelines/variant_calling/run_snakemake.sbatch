#!/usr/bin/env bash
#SBATCH --job-name=variant_workflow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=workflow-%j.log
set -euo pipefail
# Activate the workflow's software environment before submitting this file.
# Requires the Snakemake 7 legacy --cluster interface and a shared filesystem.
WORKFLOW_DIR="$(realpath "${1:?Usage: sbatch run_snakemake.sbatch WORKFLOW_DIR}")"
cd "$WORKFLOW_DIR"
mkdir -p logs
# Resolve relative configuration/sample paths without editing the Snakefile.
RUNTIME_CONFIG="$(python3 ../../prepare_paths.py snake_conf.yaml)"
snakemake --snakefile Snakefile --configfile "$RUNTIME_CONFIG" --jobs "${JOBS:-40}" \
  --cluster "sbatch --nodes=1 --ntasks=1 --cpus-per-task=${JOB_CPUS:-12} --mem=${JOB_MEM:-96G} --time=${JOB_TIME:-7-00:00:00} --output=logs/%j.log" \
  --rerun-incomplete --latency-wait 240 --restart-times 1 \
  ${SNAKEMAKE_EXTRA_ARGS:-}
