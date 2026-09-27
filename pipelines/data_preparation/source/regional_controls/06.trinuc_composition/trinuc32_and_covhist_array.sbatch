# =========================
# File: trinuc32_and_covhist_array.sbatch
# =========================
#!/usr/bin/env bash
#SBATCH --job-name=tri32_covhist
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --array=0-33%6
#SBATCH --output=log/%x.%A_%a.out
#SBATCH --error=log/%x.%A_%a.err
#SBATCH --export=ALL

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

export LC_ALL=C

# -------------------------
# Config
# -------------------------
BASE_DIR="./work/regional_controls/06.trinuc_composition"
SRC_BASE="${BASE_DIR}/baseQ13_hg38_v4/round2_chr_wise"

OUT_TRI="${BASE_DIR}/01.trinuc32"
OUT_COV="${BASE_DIR}/02.coverage_hist"

CHR_LIST=(chr{1..22} chrX)

# Make output dirs
mkdir -p \
  "${OUT_TRI}/by_chrom_tmp" \
  "${OUT_TRI}/by_chrom" \
  "${OUT_TRI}/merged" \
  "${OUT_COV}/00.by_chrom_tmp" \
  log

# -------------------------
# Methods list (auto-discover, stable order)
# -------------------------
mapfile -t METHODS < <(find "${SRC_BASE}" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort)

if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#METHODS[@]} )); then
  echo "[ERROR] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range (0..$(( ${#METHODS[@]} - 1 )))" >&2
  exit 2
fi

METHOD="${METHODS[${SLURM_ARRAY_TASK_ID}]}"
SRCDIR="${SRC_BASE}/${METHOD}"

echo "[INFO] METHOD=${METHOD}"
echo "[INFO] SRCDIR=${SRCDIR}"

# Prefer node-local tmp if available
TMPROOT="${SLURM_TMPDIR:-/tmp}/${USER}/tri32_covhist/${METHOD}_$$"
mkdir -p "${TMPROOT}"
trap 'rm -rf "${TMPROOT}"' EXIT

# -------------------------
# 32-class preset order (pyrimidine-first already in input)
# -------------------------
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

# Helper: write tri32 counts for one chr file
# Output: trinuc \t count   (always 32 lines, fixed order)
tri32_count_one() {
  local infile="$1"
  local outfile="$2"

  awk -v OFS="\t" '
    { c[$5]++ }
    END { for (k in c) print k, c[k] }
  ' "${infile}" > "${TMPROOT}/tri32.unsorted.tsv"

  awk -v OFS="\t" '
    NR==FNR { order[$1]=++i; next }
    $1 in order { val[$1]=$2 }
    END{
      # print in fixed 32 order
      for (k in order) {}
    }
  ' /dev/null /dev/null

  # Use shell to enforce order (fast enough for 32 lines)
  {
    for t in "${TRI32[@]}"; do
      v=$(awk -v key="$t" '$1==key{print $2}' "${TMPROOT}/tri32.unsorted.tsv")
      echo -e "${t}\t${v:-0}"
    done
  } > "${outfile}"
}

# Helper: write depth histogram for one chr file
# Output: depth \t count (depth>=1 only; zero handled later by you)
covhist_count_one() {
  local infile="$1"
  local outfile="$2"

  awk -v OFS="\t" '
    $4 > 0 { c[$4]++ }
    END { for (d in c) print d, c[d] }
  ' "${infile}" \
  | sort -n -k1,1 > "${outfile}"
}

# -------------------------
# Per-chrom processing
# -------------------------
for CHR in "${CHR_LIST[@]}"; do
  IN="${SRCDIR}/${CHR}.tsv"
  if [[ ! -s "${IN}" ]]; then
    echo "[WARN] Missing/empty: ${IN} (skip)"
    continue
  fi

  # trinuc32 by chrom tmp
  TRI_OUT_TMP="${OUT_TRI}/by_chrom_tmp/${METHOD}_trinuc.${CHR}.tsv"
  TRI_OUT_TMP_LOCAL="${TMPROOT}/$(basename "${TRI_OUT_TMP}")"

  tri32_count_one "${IN}" "${TRI_OUT_TMP_LOCAL}"
  mv "${TRI_OUT_TMP_LOCAL}" "${TRI_OUT_TMP}"

  # coverage hist by chrom tmp
  COV_OUT_TMP="${OUT_COV}/00.by_chrom_tmp/${METHOD}_coverage_hist.${CHR}.tsv"
  COV_OUT_TMP_LOCAL="${TMPROOT}/$(basename "${COV_OUT_TMP}")"

  covhist_count_one "${IN}" "${COV_OUT_TMP_LOCAL}"
  mv "${COV_OUT_TMP_LOCAL}" "${COV_OUT_TMP}"

  echo "[OK] ${METHOD} ${CHR}"
done

# -------------------------
# Build trinuc32 by_chrom matrix (32 rows x 23 cols)
# Output: 01.trinuc32/by_chrom/<METHOD>_trinuc.by_chrom.tsv
# -------------------------
BYCHROM_OUT="${OUT_TRI}/by_chrom/${METHOD}_trinuc.by_chrom.tsv"
BYCHROM_LOCAL="${TMPROOT}/$(basename "${BYCHROM_OUT}")"

{
  # header
  printf "trinuc"
  for CHR in "${CHR_LIST[@]}"; do
    printf "\tcount_%s" "${CHR}"
  done
  printf "\n"

  # rows
  for t in "${TRI32[@]}"; do
    printf "%s" "${t}"
    for CHR in "${CHR_LIST[@]}"; do
      f="${OUT_TRI}/by_chrom_tmp/${METHOD}_trinuc.${CHR}.tsv"
      if [[ -s "${f}" ]]; then
        v=$(awk -v key="${t}" '$1==key{print $2; exit}' "${f}")
        printf "\t%s" "${v:-0}"
      else
        printf "\t0"
      fi
    done
    printf "\n"
  done
} > "${BYCHROM_LOCAL}"

mv "${BYCHROM_LOCAL}" "${BYCHROM_OUT}"

# -------------------------
# Build trinuc32 merged (32 rows x 2 cols)
# Output: 01.trinuc32/merged/<METHOD>_trinuc.merged.tsv
# -------------------------
MERGED_OUT="${OUT_TRI}/merged/${METHOD}_trinuc.merged.tsv"
MERGED_LOCAL="${TMPROOT}/$(basename "${MERGED_OUT}")"

{
  printf "trinuc\tcount\n"
  for t in "${TRI32[@]}"; do
    sum=0
    for CHR in "${CHR_LIST[@]}"; do
      f="${OUT_TRI}/by_chrom_tmp/${METHOD}_trinuc.${CHR}.tsv"
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
} > "${MERGED_LOCAL}"

mv "${MERGED_LOCAL}" "${MERGED_OUT}"

echo "[DONE] ${METHOD}"
