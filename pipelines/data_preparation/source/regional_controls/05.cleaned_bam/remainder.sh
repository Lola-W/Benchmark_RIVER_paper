#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# FULL PATHS (all fixed here)
# -----------------------------

# Input TSV (variant list)
#./work/regional_controls/02.merged
TSV="./work/regional_controls/05.cleaned_bam/benchmark_snv_matrix.annot.with_hg19.tsv"

# Chunk VCF output dir
CHUNK_DIR="./work/regional_controls/05.cleaned_bam/chunk_100_hg38_with_down"

# Scripts (full paths)
MAKE_CHUNKS_PY="./work/regional_controls/05.cleaned_bam/make_chunk_vcfs_from_tsv.py"
WORKER_SBATCH="./work/regional_controls/05.cleaned_bam/ci_array_worker.sbatch"
AGG_SH="./work/regional_controls/05.cleaned_bam/aggregate_ci_parts.sh"

# Existing CI python (intact)
CI_PY="./work/regional_controls/05.cleaned_bam/ci_baseQ13.py"

# Reference FASTA
FASTA="./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Base tmp/out (same style as your existing)
BASE_TMP="./work/regional_controls/05.cleaned_bam/baseQ13_hg38/tmp"
BASE_OUT="./work/regional_controls/05.cleaned_bam/baseQ13_hg38/result"

# Make sure log dir exists (for worker sbatch output/error)
mkdir -p "${BASE_TMP}/_ARRAY_LOGS"

# -----------------------------
# BAMs (full paths, includes ppmseq_control)
# -----------------------------

baseDir="./work/regional_controls/05.cleaned_bam"
declare -A BAM
#BAM[nanoseq_clean]="${baseDir}/duplex_clean.NanoSeq_7614.chr1-22_X.bam"
#BAM[nanoseq_clean_6566]="${baseDir}/duplex_clean.NanoSeq_6566.chr1-22_X.bam"
#BAM[udseq_clean]="${baseDir}/duplex_clean.UDSeq_7614.chr1-22_X.bam"
#BAM[30x_udseq_clean]="${baseDir}/duplex_clean.30x_UDSeq_7614.chr1-22_X.bam"
#BAM[udseq_clean_6566]="${baseDir}/duplex_clean.UDSeq_6566.chr1-22_X.bam"
#BAM[ppmseq_clean]="${baseDir}/ppmSeq_7614.mixed_mixed.bam"
#BAM[30x_ppmseq_clean]="${baseDir}/ppmSeq_30x_7614.mixed_mixed.bam"
#BAM[100x_ppmseq_clean]="${baseDir}/ppmSeq_100x_7614.mixed_mixed.bam"
#BAM[200x_ppmseq_clean]="${baseDir}/ppmSeq_200x_7614.mixed_mixed.bam"
#BAM[ppmseq_clean_6566]="${baseDir}/ppmSeq_6566.mixed_mixed.bam"
BAM[hidef_clean]="${baseDir}/hidef_duplex.7614.7614.7614.ccs.filtered.aligned.sorted.chr1-22_X.bam"

#PLATFORMS=(nanoseq_clean nanoseq_clean_6566 udseq_clean 30x_udseq_clean udseq_clean_6566 ppmseq_clean 30x_ppmseq_clean 100x_ppmseq_clean 200x_ppmseq_clean ppmseq_clean_6566)

PLATFORMS=(hidef_clean)
# -----------------------------
# 1) Make chunk VCFs (always)
# -----------------------------
#mkdir -p "$CHUNK_DIR"
#if [[ ! -s "${CHUNK_DIR}/part_000.vcf" ]]; then
#  python3 "$MAKE_CHUNKS_PY" -i "$TSV" -o "$CHUNK_DIR" -n 100
#fi

# quick sanity: ensure 100 files exist
for i in $(seq 0 99); do
  f="${CHUNK_DIR}/part_$(printf "%03d" $i).vcf"
  if [[ ! -s "$f" ]]; then
    echo "ERROR: missing/empty chunk VCF: $f" >&2
    exit 2
  fi
done

# -----------------------------
# 2) Submit per-platform arrays + aggregation dependent job
# -----------------------------
for p in "${PLATFORMS[@]}"; do
  workdir="${BASE_TMP}/${p}"
  outdir="${BASE_OUT}/${p}/parts"
  final="${BASE_OUT}/${p}/baseQ13.big_with_downsampling.maf_ci.tsv"

  mkdir -p "$workdir" "$outdir" "$(dirname "$final")"

  echo "[*] submit array for $p"

  jid=$(sbatch \
    --job-name="ci_${p}" \
    --export=ALL,CI_FASTA="$FASTA",CI_BAM="${BAM[$p]}",CI_CI_PY="$CI_PY",CI_CHUNK_DIR="$CHUNK_DIR",CI_WORKDIR="$workdir",CI_OUTDIR="$outdir" \
    "$WORKER_SBATCH" | awk '{print $4}')

  echo "[*] array jobid for $p: $jid"
  echo "[*] submit aggregate afterok:$jid"

  sbatch \
    --job-name="ciAgg_${p}" \
     \
    --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=01:00:00 \
    --dependency="afterok:${jid}" \
    --output="${workdir}/aggregate.%j.log" --error="${workdir}/aggregate.%j.err" \
    --wrap="bash \"$AGG_SH\" \"$outdir\" \"$final\""
done

