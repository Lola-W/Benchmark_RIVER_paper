#!/usr/bin/env bash
set -euo pipefail

TEMPLATE="./work/project/code/coverage/mosdepth_template.sbatch.in"

# runner + inputs
RUNNER="./work/project/code/coverage/mosdepth_blacklist_flags.sh"
BLACKLIST="./work/project/data/filter_bam/all_blacklist.bed"
FASTA="./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# outputs
BASE_TMP="./work/project/data/coverage/tmp"
BASE_OUT="./work/project/data/coverage"

# resources
CPUS=12
MEM="96G"
TIME="3-00:00:00"

# mosdepth/run params
MAPQ=0
EXCL_FLAG=3840
THRESHOLDS="1,5,10,15,20,25,30,35,40,45,50,100,150,200,250,300,350,400,450,500,1000"

# per-platform BAMs
declare -A BAM
# BAM[illumina]="./work/project/data/illumina_hg38_bam/illumina_30x.bam"
# BAM[nanoseq]="duplex_clean.NanoSeq_7614.chr1-22_X.bam"
# BAM[udseq]="duplex_clean.30x_UDSeq_7614.chr1-22_X.bam"
# BAM[hidefseq]="hidef_duplex.7614.7614.7614.ccs.filtered.aligned.sorted.chr1-22_X.bam"
# BAM[ppmseq]="ppmSeq_30x_7614.mixed_mixed.bam"

# PLATFORMS=(illumina nanoseq udseq hidefseq ppmseq) 
# BAM[nanoseq_raw]="./inputs/alignments/7614-Cortex-L-T_Hpy_D.bam"
# BAM[udseq_raw]="./work/project/data/donwsampling/udseq/30x_7614-Cortex-L-T_D.bam"
# BAM[hidefseq_raw]="./inputs/alignments/m84137_250211_013035_s1.bam"
# BAM[ppmseq_raw]="./work/project/data/donwsampling/ppmseq/30x_7614.cram"
# PLATFORMS=(nanoseq_raw udseq_raw hidefseq_raw ppmseq_raw)

BAM[ppmseq_raw]="./work/project/data/donwsampling/ppmseq/30x_7614.cram"
PLATFORMS=(ppmseq_raw)

mkdir -p "${BASE_TMP}" "${BASE_OUT}"

for p in "${PLATFORMS[@]}"; do
  workdir="${BASE_TMP}/${p}"
  outdir="${BASE_OUT}/${p}"
  prefix="${p}.autosomesXY.noBL"
  sb="mosdepth_${p}.sbatch"

  mkdir -p "${workdir}" "${outdir}"

  sed \
    -e "s|{{PLATFORM}}|${p}|g" \
    -e "s|{{CPUS}}|${CPUS}|g" \
    -e "s|{{MEM}}|${MEM}|g" \
    -e "s|{{TIME}}|${TIME}|g" \
    -e "s|{{RUNNER}}|${RUNNER}|g" \
    -e "s|{{FASTA}}|${FASTA}|g" \
    -e "s|{{BLACKLIST}}|${BLACKLIST}|g" \
    -e "s|{{BAM}}|${BAM[$p]}|g" \
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