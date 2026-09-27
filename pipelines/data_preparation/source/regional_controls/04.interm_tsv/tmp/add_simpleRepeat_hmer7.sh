set -euo pipefail

IN="_raw_alphagenome.postqc.slim.benchmark.tsv"
OUT="raw_alphagenome.postqc.slim.benchmark.tsv"
BED_SR="simple_repeats_ucsc_hg38.sorted.bed"
BED_H7="hmers_7_and_higher.sorted.bed"

# 1) TSV -> BED (row id = NR)
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1,$2-1,$2,NR}' "$IN" \
| sort -k1,1 -k2,2n \
> _tmp.pos_with_id.bed

# 2) overlap row-ids (simple repeats)
bedtools intersect -u -a _tmp.pos_with_id.bed -b "$BED_SR" \
| cut -f4 | sort -n > _tmp.hit_simple.txt

# 3) overlap row-ids (hmer7+)
bedtools intersect -u -a _tmp.pos_with_id.bed -b "$BED_H7" \
| cut -f4 | sort -n > _tmp.hit_hmer7.txt

# 4) insert flags right after blacklist_binary, before gnomad_common
awk -v OFS="\t" '
BEGIN{
  while((getline < "_tmp.hit_simple.txt")>0) hitS[$1]=1
  while((getline < "_tmp.hit_hmer7.txt")>0) hitH[$1]=1
}
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="blacklist_binary") bbin=i
    if($i=="gnomad_common") gno=i
  }
  if(!bbin){print "ERROR: blacklist_binary not found" > "/dev/stderr"; exit 2}
  if(!gno){print "ERROR: gnomad_common not found" > "/dev/stderr"; exit 2}

  for(i=1;i<=NF;i++){
    printf "%s", $i
    if(i==bbin){
      printf "%s%s%s", OFS, "in_simple_repeats_ucsc", OFS
      printf "%s", "in_hmer7"
    }
    if(i<NF) printf "%s", OFS
    else printf "\n"
  }
  next
}
{
  s = (NR in hitS)?1:0
  h = (NR in hitH)?1:0

  for(i=1;i<=NF;i++){
    printf "%s", $i
    if(i==bbin){
      printf "%s%d%s%d", OFS, s, OFS, h
    }
    if(i<NF) printf "%s", OFS
    else printf "\n"
  }
}
' "$IN" > "$OUT"

rm -f _tmp.pos_with_id.bed _tmp.hit_simple.txt _tmp.hit_hmer7.txt

