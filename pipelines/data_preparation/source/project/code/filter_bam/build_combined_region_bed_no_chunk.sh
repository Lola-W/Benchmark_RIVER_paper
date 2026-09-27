#!/usr/bin/env bash
set -euo pipefail

BIGWIG_TO_BEDGRAPH="${BIGWIG_TO_BEDGRAPH:-bigWigToBedGraph}"
BEDTOOLS="${BEDTOOLS:-bedtools}"

OUT_BED="${1:-./work/project/data/filter_bam/all_blacklist.bed}"
OUT_DIR="$(dirname "$OUT_BED")"
mkdir -p "$OUT_DIR"

if [[ ! -x "$BIGWIG_TO_BEDGRAPH" ]]; then
  echo "ERROR: bigWigToBedGraph not executable: $BIGWIG_TO_BEDGRAPH" >&2
  exit 1
fi

if [[ ! -x "$BEDTOOLS" ]]; then
  echo "ERROR: bedtools not executable: $BEDTOOLS" >&2
  exit 1
fi

tmpdir="$(mktemp -d "$OUT_DIR/tmp.combine_regions.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

BW_FILES=(
  "./resources/hg38/bw_mask/hg38.segdup.bw"
  "./resources/hg38/bw_mask/hg38_UCSC_gaps.bw"
  "./resources/hg38/bw_mask/satellites_rDNA.bw"
  "./resources/hg38/k50.Umap.MultiTrackMappability.bw"
  "./resources/hg38/bw_mask/ENCODE_hg38-blacklist.v2.bw"
)

BW_THRESHOLDS=(
  "gte0.1"
  "gte0.1"
  "gte0.1"
  "lt0.4"
  "gte0.1"
)

EXTRA_BEDS=(
  "./work/regional_controls/simple_repeats_ucsc_hg38.sorted.bed"
  "./work/regional_controls/hmers_7_and_higher.sorted.bed"
)

filter_bigwig_to_bed() {
  local bw_file="$1"
  local threshold="$2"
  local out_file="$3"
  local op value

  if [[ "$threshold" =~ ^(lt|lte|gt|gte)([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)$ ]]; then
    op="${BASH_REMATCH[1]}"
    value="${BASH_REMATCH[2]}"
  else
    echo "ERROR: unsupported threshold format: $threshold" >&2
    exit 1
  fi

  "$BIGWIG_TO_BEDGRAPH" "$bw_file" /dev/stdout \
    | awk -v OFS='\t' -v op="$op" -v value="$value" '
        function keep(x) {
          if (op == "lt")  return (x < value)
          if (op == "lte") return (x <= value)
          if (op == "gt")  return (x > value)
          if (op == "gte") return (x >= value)
          return 0
        }
        keep($4) { print $1, $2, $3 }
      ' > "$out_file"
}

echo "Writing temporary region lists into: $tmpdir"

for i in "${!BW_FILES[@]}"; do
  bw="${BW_FILES[$i]}"
  threshold="${BW_THRESHOLDS[$i]}"
  out="$tmpdir/bw_${i}.bed"

  if [[ ! -f "$bw" ]]; then
    echo "ERROR: missing BigWig: $bw" >&2
    exit 1
  fi

  echo "Filtering $(basename "$bw") with threshold $threshold"
  filter_bigwig_to_bed "$bw" "$threshold" "$out"
done

for i in "${!EXTRA_BEDS[@]}"; do
  bed="${EXTRA_BEDS[$i]}"
  out="$tmpdir/extra_${i}.bed"

  if [[ ! -f "$bed" ]]; then
    echo "ERROR: missing BED: $bed" >&2
    exit 1
  fi

  echo "Adding BED $(basename "$bed")"
  awk 'BEGIN { OFS="\t" } !/^#/ && NF >= 3 { print $1, $2, $3 }' "$bed" > "$out"
done

cat "$tmpdir"/*.bed \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | "$BEDTOOLS" merge -i - > "$OUT_BED"

echo "Done."
echo "Output BED: $OUT_BED"
