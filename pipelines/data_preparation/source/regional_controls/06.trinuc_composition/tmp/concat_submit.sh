BASE="baseQ13_hg38_v4/round2_chr_wise"

for d in ${BASE}/*; do
  [ -d "$d" ] || continue
  method=$(basename "$d")

  sbatch --job-name="concat_${method}" \
         --export=METHOD="$method",SRCDIR="$d" \
         concat_ibp.sbatch
done
