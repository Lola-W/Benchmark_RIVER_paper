# 1) list files: ref first, then others sorted
ref="ref_hg38_trinuc.merged.tsv"
others=$(ls *_trinuc.merged.tsv | grep -v "^${ref}$" | sort)

# 2) header
{
  printf "trinuc\tref_hg38"
  for f in $others; do
    printf "\t%s" "${f%_trinuc.merged.tsv}"
  done
  printf "\n"
} > trinuc32_all_methods_wide.tsv

# 3) body (paste trinuc + count cols)
paste \
  <(awk 'NR>1{print $1}' "$ref") \
  <(awk 'NR>1{print $2}' "$ref") \
  $(for f in $others; do echo "<(awk 'NR>1{print \$2}' $f)"; done) \
  >> trinuc32_all_methods_wide.tsv
