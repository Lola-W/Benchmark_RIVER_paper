#!/usr/bin/env bash
#SBATCH --job-name=hidefseq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=hidefseq-%j.log
set -euo pipefail
# Submit from this directory after editing params.yaml and nextflow.config.
# HiDEF-seq release 3.1; commit recorded by the retained study run.
RUNTIME_PARAMS="$(python3 ../prepare_paths.py params.yaml)"
nextflow -config nextflow.config run evronylab/hidef-seq \
  -r 15c522072cbc9ea04e521c75ac047029c53c2d23 \
  -params-file "$RUNTIME_PARAMS" -resume -work-dir "${WORK_DIR:-./work}" -with-report
