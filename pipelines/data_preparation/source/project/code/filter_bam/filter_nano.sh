set -euo pipefail

BAM="./inputs/alignments/7614-Cortex-L-T_Hpy_D.bam"
SCRIPT="./work/project/code/filter_bam/filter_duplex_umi_families.py"
OUT="./work/project/data/filter_bam/nano"
mkdir -p "$OUT/chrom"

samtools idxstats "$BAM" | cut -f1 | grep -v '^\*$' > "$OUT/chrom.list"

cat "$OUT/chrom.list" | xargs -I{} -P8 bash -lc '
  samtools view -b "'"$BAM"'" "{}" -o "'"$OUT"'"/chrom/{}.bam
  python "'"$SCRIPT"'" "'"$OUT"'"/chrom/{}.bam "'"$OUT"'"/chrom/{}.umi_counts.tsv "'"$OUT"'"/chrom/{}.cleaned.bam
'

samtools merge -@8 -f "$OUT/cleaned.merged.bam" "$OUT"/chrom/*.cleaned.bam
