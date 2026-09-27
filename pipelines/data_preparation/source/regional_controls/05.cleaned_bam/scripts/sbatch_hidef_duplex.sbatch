#!/usr/bin/env bash
#SBATCH --job-name=hidef_duplex
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --export=ALL
#SBATCH --output=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/hidef_duplex.%j.log
#SBATCH --error=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/hidef_duplex.%j.err

set -euo pipefail

# -----------------------------
# Env
# -----------------------------
set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

# -----------------------------
# Fixed paths (all variables in-file as requested)
# -----------------------------
OUT_BASEDIR="./work/regional_controls/05.cleaned_bam"
IN_BAM="./inputs/alignments/7614.7614.7614.ccs.filtered.aligned.sorted.bam"

EXCL_FLAGS=3840
CHRS=(chr{1..22} chrX)

mkdir -p "${OUT_BASEDIR}/tmp/_CHR_LOGS"
mkdir -p "${OUT_BASEDIR}/hidef_clean"
mkdir -p "${OUT_BASEDIR}/tmp/hidef_clean_work"

BASENAME="$(basename "${IN_BAM}")"
STEM="${BASENAME%.cram}"
STEM="${STEM%.bam}"

WORKDIR="${OUT_BASEDIR}/tmp/hidef_clean_work/${STEM}"
mkdir -p "${WORKDIR}"

ZMW_LIST="${WORKDIR}/duplex_zmw.${STEM}.chr1-22_X.txt"
OUT_UNSORTED="${WORKDIR}/hidef_duplex.${STEM}.chr1-22_X.unsorted.bam"
OUT_BAM="${OUT_BASEDIR}/hidef_clean/hidef_duplex.${STEM}.chr1-22_X.bam"

echo "[INFO] IN_BAM=${IN_BAM}"
echo "[INFO] OUT_BAM=${OUT_BAM}"
echo "[INFO] EXCL_FLAGS=${EXCL_FLAGS}"
echo "[INFO] CHRS=${CHRS[*]}"

# -----------------------------
# 1) Build duplex ZMW list across chr1-22,chrX after -F 3840
#    duplex criterion: >=1 fwd AND >=1 rev over the aggregated stream
#    qname assumed like: movie/zmw/ccs/fwd (strand label in 4th token)
# -----------------------------
samtools view -@ "${SLURM_CPUS_PER_TASK}" -F "${EXCL_FLAGS}" "${IN_BAM}" "${CHRS[@]}" \
| awk '
{
  n=split($1,a,"/")
  if(n<2) next
  zm=a[2]
  q4=(n>=4?a[4]:"NA")
  if(q4=="fwd") f[zm]++
  else if(q4=="rev") r[zm]++
}
END{
  for(zm in f){
    if((f[zm]+0)>=1 && (r[zm]+0)>=1) print zm
  }
  for(zm in r){
    if((f[zm]+0)>=1 && (r[zm]+0)>=1) print zm
  }
}' | sort -u > "${ZMW_LIST}"

echo "[INFO] duplex_zmw_count=$(wc -l < "${ZMW_LIST}")"

# -----------------------------
# 2) Emit only reads whose ZMW is duplex; write BAM (unsorted)
# -----------------------------
samtools view -h -@ "${SLURM_CPUS_PER_TASK}" -F "${EXCL_FLAGS}" "${IN_BAM}" "${CHRS[@]}" \
| awk -v L="${ZMW_LIST}" '
BEGIN{
  while((getline < L)>0) keep[$1]=1
  close(L)
}
/^@/ { print; next }
{
  n=split($1,a,"/")
  if(n<2) next
  zm=a[2]
  if(keep[zm]) print
}
' \
| samtools view -@ "${SLURM_CPUS_PER_TASK}" -b -o "${OUT_UNSORTED}" -

# -----------------------------
# 3) Sort + index (final output under hidef_clean/)
# -----------------------------
TMP_SORT="${WORKDIR}/sorted.bam"
samtools sort \
  -@ "${SLURM_CPUS_PER_TASK}" \
  -m 6G \
  -o "${TMP_SORT}" \
  "${OUT_UNSORTED}"

mv "${TMP_SORT}" "${OUT_BAM}"
samtools index -@ "${SLURM_CPUS_PER_TASK}" "${OUT_BAM}"

rm -f "${OUT_UNSORTED}"

echo "[DONE] ${OUT_BAM}"
echo "[DONE] ${OUT_BAM}.bai"
