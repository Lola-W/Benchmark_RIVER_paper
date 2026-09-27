#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# FULL PATHS (all fixed here)
# -----------------------------

# Input TSV (variant list)
#./work/regional_controls/02.merged
TSV="./work/regional_controls/02.merged/benchmark_snv_matrix.annot.with_hg19.tsv"

# Chunk VCF output dir
CHUNK_DIR="./work/regional_controls/02.merged/chunk_100_hg38_with_down"

# Scripts (full paths)
MAKE_CHUNKS_PY="./work/regional_controls/02.merged/make_chunk_vcfs_from_tsv.py"
WORKER_SBATCH="./work/regional_controls/02.merged/ci_array_worker.sbatch"
AGG_SH="./work/regional_controls/02.merged/aggregate_ci_parts.sh"

# Existing CI python (intact)
CI_PY="./work/regional_controls/02.merged/ci_baseQ13.py"

# Reference FASTA
FASTA="./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# Base tmp/out (same style as your existing)
BASE_TMP="./work/regional_controls/02.merged/down_baseQ13_hg38/tmp"
BASE_OUT="./work/regional_controls/02.merged/down_baseQ13_hg38/result"

# Make sure log dir exists (for worker sbatch output/error)
mkdir -p "${BASE_TMP}/_ARRAY_LOGS"

# -----------------------------
# BAMs (full paths, includes ppmseq_control)
# -----------------------------
declare -A BAM
#BAM[illumina]="./work/project/data/illumina_hg38_bam/7614_illumina_brain_somatic.bam"
#BAM[nanoseq]="./inputs/alignments/7614-Cortex-L-T_Hpy_D.bam"
#BAM[nanoseq_6566]="./inputs/alignments/S6566_Hpy.bam"
#BAM[udseq]="./inputs/alignments/7614-Cortex-L-T_D.bam"
#BAM[udseq_6566]="./inputs/alignments/6566.bam"
#BAM[hidefseq]="./inputs/alignments/m84137_250211_013035_s1.bam"
#BAM[ppmseq]="./inputs/alignments/422077-25-5853-DNA-1-ppm0022-CACAACATATCAGAT.cram"
#BAM[ppmseq_6566]="./inputs/alignments/422073-25-5852-DNA-1-ppm0021-CGCATCCTCACAGAT.cram"

BAM[ppmseq_30x]="./work/project/data/donwsampling/ppmseq/30x_7614.cram"
BAM[ppmseq_100x]="./work/project/data/donwsampling/ppmseq/100x_7614.cram"
BAM[ppmseq_200x]="./work/project/data/donwsampling/ppmseq/200x_7614.cram"

BAM[illumina_30x]="./work/project/data/illumina_hg38_bam/illumina_30x.bam"
BAM[illumina_100x]="./work/project/data/illumina_hg38_bam/illumina_100x.bam"
BAM[illumina_200x]="./work/project/data/illumina_hg38_bam/illumina_200x.bam"

BAM[udseq_30x]="./work/project/data/donwsampling/udseq/30x_7614-Cortex-L-T_D.bam"



# -----------------------------
# PTA single-neuron BAMs
# -----------------------------
#PTA_BAM_DIR="./inputs/pta_hg38"

#for i in {1..16}; do
#  key="PTA_single_neuron_$(printf "%02d" $i)"
#  BAM["$key"]="${PTA_BAM_DIR}/7614_single_neuron_${i}.bam"
#done

#PLATFORMS=(illumina nanoseq nanoseq_6566 udseq udseq_6566 hidefseq ppmseq ppmseq_6566)
PLATFORMS=(ppmseq_30x ppmseq_100x ppmseq_200x illumina_30x illumina_100x illumina_200x udseq_30x)
#PLATFORMS+=( $(printf "PTA_single_neuron_%02d " {1..16}) )

# -----------------------------
# 1) Make chunk VCFs (always)
# -----------------------------
#mkdir -p "$CHUNK_DIR"
#if [[ ! -s "${CHUNK_DIR}/part_000.vcf" ]]; then
#  python3 "$MAKE_CHUNKS_PY" -i "$TSV" -o "$CHUNK_DIR" -n 100
#fi

# quick sanity: ensure 100 files exist
#for i in $(seq 0 99); do
#  f="${CHUNK_DIR}/part_$(printf "%03d" $i).vcf"
#  if [[ ! -s "$f" ]]; then
#    echo "ERROR: missing/empty chunk VCF: $f" >&2
#    exit 2
#  fi
#done

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

