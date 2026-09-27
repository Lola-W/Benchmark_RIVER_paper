#!/usr/bin/env bash
set -euo pipefail

TEMPLATE="./work/project/code/coverage/mosdepth_template.sbatch.in"

# runner + inputs (reuse bulk)
RUNNER="./work/project/code/coverage/mosdepth_blacklist_flags.sh"
BLACKLIST="./work/project/data/filter_bam/all_blacklist.bed"
FASTA="./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# outputs (reuse bulk base)
BASE_TMP="./work/project/data/coverage/tmp"
BASE_OUT="./work/project/data/coverage"

# resources (reuse bulk)
CPUS=12
MEM="96G"
TIME="3-00:00:00"

# mosdepth/run params (reuse bulk)
MAPQ=0
EXCL_FLAG=3840
THRESHOLDS="1,5,10,15,20,25,30,35,40,45,50,100,150,200,250,300,350,400,450,500,1000"

# PTA BAM dir (EDIT)
PTA_BAM_DIR="./inputs/pta_hg38"

mkdir -p "${BASE_TMP}" "${BASE_OUT}"

PTA_TMP="${BASE_TMP}/pta"
PTA_OUT="${BASE_OUT}/pta"
mkdir -p "${PTA_TMP}" "${PTA_OUT}"

find "${PTA_BAM_DIR}" -maxdepth 1 -type f -name "*.bam" | sort | while read -r bam; do
  sample="$(basename "$bam" .bam)"

  workdir="${PTA_TMP}/${sample}"
  outdir="${PTA_OUT}"
  prefix="${sample}.autosomesXY.noBL"
  sb="${workdir}/mosdepth_${sample}.sbatch"

  mkdir -p "${workdir}" "${outdir}"

  sed \
    -e "s|{{PLATFORM}}|pta_${sample}|g" \
    -e "s|{{CPUS}}|${CPUS}|g" \
    -e "s|{{MEM}}|${MEM}|g" \
    -e "s|{{TIME}}|${TIME}|g" \
    -e "s|{{RUNNER}}|${RUNNER}|g" \
    -e "s|{{FASTA}}|${FASTA}|g" \
    -e "s|{{BLACKLIST}}|${BLACKLIST}|g" \
    -e "s|{{BAM}}|${bam}|g" \
    -e "s|{{WORKDIR}}|${workdir}|g" \
    -e "s|{{OUTDIR}}|${outdir}|g" \
    -e "s|{{PREFIX}}|${prefix}|g" \
    -e "s|{{MAPQ}}|${MAPQ}|g" \
    -e "s|{{EXCL_FLAG}}|${EXCL_FLAG}|g" \
    -e "s|{{THRESHOLDS}}|${THRESHOLDS}|g" \
    "${TEMPLATE}" > "${sb}"

  echo "[*] Submitting ${sb}"
  sbatch "${sb}"
done
