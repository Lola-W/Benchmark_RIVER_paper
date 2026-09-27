#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# hg19 CI runner (independent from hg38)
# ============================================================

# -----------------------------
# INPUT / REFERENCES
# -----------------------------

# TSV with placeholders (row order preserved)
#./work/regional_controls/02.merged
TSV="./work/regional_controls/02.merged/benchmark_snv_matrix.annot.with_hg19.tsv"

# hg19 reference fasta
FASTA="./work/regional_controls/human_g1k_v37_decoy.fasta"

# -----------------------------
# CHUNKING
# -----------------------------
CHUNK_DIR="./work/regional_controls/02.merged/chunk_100_hg19_with_down"
MAKE_CHUNKS_PY="./work/regional_controls/02.merged/make_chunk_vcfs_from_tsv_hg19.py"

# -----------------------------
# CI PIPELINE SCRIPTS
# -----------------------------
CI_PY="./work/regional_controls/02.merged/ci_baseQ13.py"
WORKER_SBATCH="./work/regional_controls/02.merged/ci_array_worker_hg19.sbatch"
AGG_SH="./work/regional_controls/02.merged/aggregate_ci_parts.sh"

# -----------------------------
# OUTPUT BASE
# -----------------------------
BASE_TMP="./work/regional_controls/02.merged/baseQ13_hg19/tmp"
BASE_OUT="./work/regional_controls/02.merged/baseQ13_hg19/result"
LOG_DIR="${BASE_TMP}/_ARRAY_LOGS"

mkdir -p "$LOG_DIR"

# -----------------------------
# BAMs (hg19)
# naming rule enforced in keys
# -----------------------------
declare -A BAM

BAM[Illumina_7669]="./inputs/alignments/7669-Br-L-P-1.bam"

# ---- PTA control (hg19) ----
BAM[PTA_6566_1]="./inputs/alignments/Young-6566-2sc3.merged.bam"
BAM[PTA_6566_2]="./inputs/alignments/Young-6566-2sc22.merged.bam"
BAM[PTA_8420_1]="./inputs/alignments/Old-8420-SC06.merged.bam"
BAM[PTA_8420_2]="./inputs/alignments/Old-8420-SC15.merged.bam"

# ---- Bulk Illumina (hg19, Lrg) ----
BAM[Illumina_B-L-F-Lrg]="./inputs/alignments/7614-B-L-F-3-Lrg.bam"
BAM[Illumina_B-L-O-Lrg]="./inputs/alignments/7614-B-L-O-5-Lrg.bam"
BAM[Illumina_B-L-P-Lrg]="./inputs/alignments/7614-B-L-P-4-Lrg.bam"
BAM[Illumina_B-L-PF-Lrg]="./inputs/alignments/7614-B-L-PF-2-Lrg.bam"
BAM[Illumina_B-L-T-Lrg]="./inputs/alignments/7614-B-L-T-1-Lrg.bam"

BAM[Illumina_B-R-F-Lrg]="./inputs/alignments/7614-B-R-F-3-Lrg.bam"
BAM[Illumina_B-R-O-Lrg]="./inputs/alignments/7614-B-R-O-5-Lrg.bam"
BAM[Illumina_B-R-P-Lrg]="./inputs/alignments/7614-B-R-P-4-Lrg.bam"
BAM[Illumina_B-R-PF-Lrg]="./inputs/alignments/7614-B-R-PF-2-Lrg.bam"
BAM[Illumina_B-R-T-Lrg]="./inputs/alignments/7614-B-R-T-1-Lrg.bam"

# ============================================================
# 1) Make hg19 chunk VCFs (once)
# ============================================================
mkdir -p "$CHUNK_DIR"

if [[ ! -s "${CHUNK_DIR}/part_000.vcf" ]]; then
  echo "[*] making hg19 chunk VCFs"
  python3 "$MAKE_CHUNKS_PY" \
    -i "$TSV" \
    -o "$CHUNK_DIR" \
    -n 100
fi

# sanity
for i in $(seq 0 99); do
  f="${CHUNK_DIR}/part_$(printf "%03d" $i).vcf"
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] missing chunk: $f" >&2
    exit 2
  fi
done

# ============================================================
# 2) Submit per-BAM array jobs + aggregation
# ============================================================
for key in "${!BAM[@]}"; do
  bam="${BAM[$key]}"
  workdir="${BASE_TMP}/${key}"
  outdir="${BASE_OUT}/${key}/parts"
  final="${BASE_OUT}/${key}/baseQ13.big_with_downsampling.hg19.maf_ci.tsv"

  mkdir -p "$workdir" "$outdir" "$(dirname "$final")"

  echo "[*] submit hg19 array: $key"

  jid=$(sbatch \
    --job-name="ci_hg19_${key}" \
    --export=ALL,CI_FASTA="$FASTA",CI_BAM="$bam",CI_CI_PY="$CI_PY",CI_CHUNK_DIR="$CHUNK_DIR",CI_WORKDIR="$workdir",CI_OUTDIR="$outdir" \
    "$WORKER_SBATCH" | awk '{print $4}')

  echo "    jobid=${jid}"

  sbatch \
    --job-name="ciAgg_hg19_${key}" \
     \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=01:00:00 \
    --dependency="afterok:${jid}" \
    --output="${workdir}/aggregate.%j.log" \
    --error="${workdir}/aggregate.%j.err" \
    --wrap="bash \"$AGG_SH\" \"$outdir\" \"$final\""
done

echo "[OK] hg19 submission done"

