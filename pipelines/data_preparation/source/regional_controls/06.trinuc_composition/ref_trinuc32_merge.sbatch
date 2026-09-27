#!/usr/bin/env bash
#SBATCH --job-name=ref_tri32_merge
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=log/%x.%j.out
#SBATCH --error=log/%x.%j.err
#SBATCH --export=ALL

set -euo pipefail

export LC_ALL=C

BASE="./work/regional_controls/06.trinuc_composition"

TMPDIR="${BASE}/01.trinuc32/by_chrom_tmp"
BYDIR="${BASE}/01.trinuc32/by_chrom"
MDIR="${BASE}/01.trinuc32/merged"

mkdir -p "${BYDIR}" "${MDIR}" log

CHR_LIST=(chr{1..22} chrX)
TRI32=(
  ACA ACC ACG ACT
  CCA CCC CCG CCT
  GCA GCC GCG GCT
  TCA TCC TCG TCT
  ATA ATC ATG ATT
  CTA CTC CTG CTT
  GTA GTC GTG GTT
  TTA TTC TTG TTT
)

# 1) by_chrom matrix
OUT_BY="${BYDIR}/ref_hg38_trinuc.by_chrom.tsv"
OUT_BY_TMP="${OUT_BY}.tmp"

{
  printf "trinuc"
  for c in "${CHR_LIST[@]}"; do
    printf "\tcount_%s" "${c}"
  done
  printf "\n"

  for t in "${TRI32[@]}"; do
    printf "%s" "${t}"
    for c in "${CHR_LIST[@]}"; do
      f="${TMPDIR}/ref_hg38_trinuc.${c}.tsv"
      if [[ ! -s "${f}" ]]; then
        printf "\t0"
      else
        v=$(awk -v key="${t}" '$1==key{print $2; exit}' "${f}")
        printf "\t%s" "${v:-0}"
      fi
    done
    printf "\n"
  done
} > "${OUT_BY_TMP}"

mv "${OUT_BY_TMP}" "${OUT_BY}"
echo "[DONE] by_chrom -> ${OUT_BY}"

# 2) merged (sum across chr)
OUT_M="${MDIR}/ref_hg38_trinuc.merged.tsv"
OUT_M_TMP="${OUT_M}.tmp"

{
  printf "trinuc\tcount\n"
  for t in "${TRI32[@]}"; do
    sum=0
    for c in "${CHR_LIST[@]}"; do
      f="${TMPDIR}/ref_hg38_trinuc.${c}.tsv"
      if [[ -s "${f}" ]]; then
        v=$(awk -v key="${t}" '$1==key{print $2; exit}' "${f}")
        v=${v:-0}
      else
        v=0
      fi
      sum=$((sum + v))
    done
    printf "%s\t%d\n" "${t}" "${sum}"
  done
} > "${OUT_M_TMP}"

mv "${OUT_M_TMP}" "${OUT_M}"
echo "[DONE] merged -> ${OUT_M}"
