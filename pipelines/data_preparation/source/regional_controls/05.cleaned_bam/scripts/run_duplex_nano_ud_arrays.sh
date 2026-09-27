#!/usr/bin/env bash
set -euo pipefail

# =========================
# CONFIG
# =========================
OUT_BASEDIR=./work/regional_controls/05.cleaned_bam
PY_SCRIPT=./work/regional_controls/05.cleaned_bam/scripts/duplex_nano_ud_clean_by_chrom.py
SBATCH_SCRIPT=./work/regional_controls/05.cleaned_bam/scripts/sbatch_duplex_nano_ud_clean_by_chrom.sbatch

# =========================
# BAM LIST (EDIT THIS)
# =========================

BAMS=(
  "./inputs/alignments/7614-Cortex-L-T_Hpy_D.bam"
  "./inputs/alignments/S6566_Hpy.bam"
  "./inputs/alignments/7614-Cortex-L-T_D.bam"
  "./work/project/data/donwsampling/udseq/30x_7614-Cortex-L-T_D.bam"
  "./inputs/alignments/6566.bam"
)

mkdir -p "${OUT_BASEDIR}"

# =========================
# SUBMIT
# =========================
for IN_BAM in "${BAMS[@]}"; do
  if [[ ! -s "${IN_BAM}" ]]; then
    echo "ERROR: missing BAM: ${IN_BAM}" >&2
    exit 1
  fi
  if [[ ! -s "${IN_BAM}.bai" ]] && [[ ! -s "${IN_BAM%.bam}.bai" ]]; then
    echo "ERROR: missing BAM index for: ${IN_BAM}" >&2
    exit 1
  fi

  export IN_BAM
  export OUT_BASEDIR
  export PY_SCRIPT

  JOBID=$(sbatch "${SBATCH_SCRIPT}" | awk '{print $4}')
  echo "Submitted: ${JOBID}  BAM=$(basename "${IN_BAM}")"
done
