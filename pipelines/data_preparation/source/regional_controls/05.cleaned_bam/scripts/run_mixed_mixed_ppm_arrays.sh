#!/usr/bin/env bash
set -euo pipefail

OUT_BASEDIR="./work/regional_controls/05.cleaned_bam/ppmSeq_mixedmixed"
SBATCH_SCRIPT="./work/regional_controls/05.cleaned_bam/scripts/sbatch_mixed_mixed_ppm_clean_by_chrom.sbatch"

declare -A BAM

BAM[ppmseq_30x]="./work/project/data/donwsampling/ppmseq/30x_7614.cram"
BAM[ppmseq_100x]="./work/project/data/donwsampling/ppmseq/100x_7614.cram"
BAM[ppmseq_200x]="./work/project/data/donwsampling/ppmseq/200x_7614.cram"
BAM[ppmseq]="./inputs/alignments/422077-25-5853-DNA-1-ppm0022-CACAACATATCAGAT.cram"
BAM[ppmseq_6566]="./inputs/alignments/422073-25-5852-DNA-1-ppm0021-CGCATCCTCACAGAT.cram"

mkdir -p "${OUT_BASEDIR}"

for key in "${!BAM[@]}"; do
  IN_CRAM="${BAM[$key]}"

  if [[ ! -s "${IN_CRAM}" ]]; then
    echo "ERROR: missing CRAM: ${IN_CRAM}"
    exit 1
  fi

  if [[ ! -s "${IN_CRAM}.crai" ]]; then
    echo "ERROR: missing CRAM index: ${IN_CRAM}.crai"
    exit 1
  fi

  export IN_CRAM
  export OUT_BASEDIR

  JOBID=$(sbatch "${SBATCH_SCRIPT}" | awk '{print $4}')
  echo "Submitted ${JOBID} → ${key}"
done
