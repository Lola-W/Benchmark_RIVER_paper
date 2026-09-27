#!/usr/bin/env bash
set -euo pipefail

# ---- fixed paths (copied from your runner) ----
CHUNK_DIR="./work/regional_controls/05.cleaned_bam/chunk_100_hg38_with_down"
WORKER_SBATCH="./work/regional_controls/05.cleaned_bam/ci_array_worker.sbatch"
AGG_SH="./work/regional_controls/05.cleaned_bam/aggregate_ci_parts.sh"
CI_PY="./work/regional_controls/05.cleaned_bam/ci_baseQ13.py"
FASTA="./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"

baseDir="./work/regional_controls/05.cleaned_bam"
BAM_nanoseq_clean="${baseDir}/duplex_clean.NanoSeq_7614.chr1-22_X.bam"

BASE_TMP="${baseDir}/baseQ13_hg38/tmp"
BASE_OUT="${baseDir}/baseQ13_hg38/result"

p="nanoseq_clean"
workdir="${BASE_TMP}/${p}"
outdir="${BASE_OUT}/${p}/parts"
final="${BASE_OUT}/${p}/baseQ13.big_with_downsampling.maf_ci.tsv"

mkdir -p "$workdir" "$outdir" "$(dirname "$final")" "${workdir}/_ARRAY_LOGS"

# ---- sanity: chunk vcf 존재 확인 (재생성 안함) ----
for i in $(seq 0 99); do
  f="${CHUNK_DIR}/part_$(printf "%03d" $i).vcf"
  [[ -s "$f" ]] || { echo "ERROR: missing/empty $f" >&2; exit 2; }
done

# ---- (optional) wipe only the failed outputs so aggregate doesn't pick stale 0-byte ----
rm -f \
  "${outdir}/part_066.tsv" \
  "${outdir}/part_067.tsv" \
  "${outdir}/part_068.tsv" \
  "${outdir}/part_069.tsv" \
  "${outdir}/part_070.tsv"

echo "[*] resubmit nanoseq_clean array 66-70 only"

jid=$(sbatch \
  --job-name="ci_${p}_rerun_66_70" \
  --array=66-70 \
  --export=ALL,CI_FASTA="$FASTA",CI_BAM="$BAM_nanoseq_clean",CI_CI_PY="$CI_PY",CI_CHUNK_DIR="$CHUNK_DIR",CI_WORKDIR="$workdir",CI_OUTDIR="$outdir" \
  "$WORKER_SBATCH" | awk '{print $4}')

echo "[*] array jobid: $jid"
echo "[*] submit aggregate afterok:$jid"

sbatch \
  --job-name="ciAgg_${p}_rerun_66_70" \
   \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G --time=01:00:00 \
  --dependency="afterok:${jid}" \
  --output="${workdir}/aggregate.rerun66_70.%j.log" --error="${workdir}/aggregate.rerun66_70.%j.err" \
  --wrap="bash \"$AGG_SH\" \"$outdir\" \"$final\""
