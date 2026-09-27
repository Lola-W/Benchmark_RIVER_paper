#!/usr/bin/env bash
set -euo pipefail

# reorder_all_blocks.sh
# - Reorders the TSV into:
#   1-4: chr,pos,ref,alt
#   Callset block (with call_ prefix; ppmSeq call order fixed; PTA summary inserted)
#   Blacklist/annotation block
#   hg19 block
#   BAM quant block (category priority: benchmark > HG19_ILLUMINA_B_*_LRG > negative controls; then method order)
#   Germline-like flags (optionally drop the final any-bulk flag)
#
# Usage:
#   bash reorder_all_blocks.sh -i IN.tsv -o OUT.tsv
#   bash reorder_all_blocks.sh -i IN.tsv -o OUT.tsv --keep_any_bulk

usage() {
  cat >&2 <<'EOF'
Usage:
  bash reorder_all_blocks.sh -i INPUT.tsv -o OUTPUT.tsv [--keep_any_bulk]

Notes:
  - Default drops "is_germline_any_bulk_upper_ci" (your col 359).
EOF
  exit 1
}

in_tsv=""
out_tsv=""
keep_any_bulk=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) in_tsv="$2"; shift 2 ;;
    -o) out_tsv="$2"; shift 2 ;;
    --keep_any_bulk) keep_any_bulk=1; shift 1 ;;
    -h|--help) usage ;;
    *) echo "[ERROR] unknown arg: $1" >&2; usage ;;
  esac
done

[[ -z "$in_tsv" || -z "$out_tsv" ]] && usage
[[ ! -f "$in_tsv" ]] && { echo "[ERROR] not found: $in_tsv" >&2; exit 2; }

gawk -v FS='\t' -v OFS='\t' -v KEEP_ANY_BULK="$keep_any_bulk" '
function die(msg){ print "[ERROR] " msg > "/dev/stderr"; exit 2 }

function push(arr_name, v,    n){
  n = ++arrN[arr_name]
  arr[arr_name, n] = v
}

function add_callset(v){ push("CALLSET", v) }
function add_black(v){  push("BLACK", v) }
function add_hg19(v){   push("HG19", v) }
function add_flag(v){   push("FLAGS", v) }

function is_bam_quant(col){
  return (col ~ /_(REF_COUNT|ALT_COUNT|MAF|LOWER_CI|UPPER_CI|IS_GREATER)$/)
}

function bam_category(col){
  # 1) same sample benchmark (default)
  # 2) same sample different region: HG19_ILLUMINA_B_*_LRG
  # 3) different sample negative control: contains 6566/7669/8420
  if(col ~ /^HG19_ILLUMINA_B_.*_LRG_/ ) return 2
  if(col ~ /(6566|7669|8420)/) return 3
  return 1
}

function bam_method(col){
  # method priority: ILLUMINA, PPMSEQ, UDSEQ, NANOSEQ, HIDEFSEQ, PTA, AmpliSeq
  if(col ~ /^ILLUMINA_/) return "ILLUMINA"
  if(col ~ /^PPMSEQ_/) return "PPMSEQ"
  if(col ~ /^UDSEQ_/) return "UDSEQ"
  if(col ~ /^NANOSEQ_/) return "NANOSEQ"
  if(col ~ /^HIDEFSEQ_/) return "HIDEFSEQ"
  if(col ~ /^PTA_SINGLE_NEURON_/) return "PTA"

  # hg19-lifted controls still map to underlying method
  if(col ~ /^HG19_ILLUMINA_/) return "ILLUMINA"
  if(col ~ /^HG19_PTA_/) return "PTA"

  return "ZZZ"
}

function method_rank(m){
  if(m=="ILLUMINA") return 1
  if(m=="PPMSEQ")   return 2
  if(m=="UDSEQ")    return 3
  if(m=="NANOSEQ")  return 4
  if(m=="HIDEFSEQ") return 5
  if(m=="PTA")      return 6
  if(m=="AMPLISEQ") return 7
  return 99
}

function bam_down_rank(col){
  # within a method: original, 30x, 100x, 200x (missing -> just those present)
  if(col ~ /_200x_/) return 4
  if(col ~ /_100x_/) return 3
  if(col ~ /_30x_/)  return 2
  return 1
}

function bam_pta_rank(col,    m){
  # PTA_SINGLE_NEURON_01..16 ordering
  if(match(col, /^PTA_SINGLE_NEURON_([0-9]{2})_/, m)) return (m[1] + 0)
  # HG19_PTA_6566_1 / HG19_PTA_6566_2 / HG19_PTA_8420_1 / _2 : keep numeric suffix
  if(match(col, /^HG19_PTA_[0-9]+_([0-9]+)_/, m)) return (m[1] + 0)
  return 999
}

function bam_sort_key(col,    cat, m, mr, dr, pr){
  cat = bam_category(col)
  m   = bam_method(col)
  mr  = method_rank(m)
  dr  = bam_down_rank(col)
  pr  = bam_pta_rank(col)

  # Key order:
  #   category (1,2,3) -> method rank -> PTA idx (if any) -> downsample rank -> lexical (stable)
  return sprintf("%d|%02d|%03d|%02d|%s", cat, mr, pr, dr, col)
}

BEGIN{
  # --- callset order you want (method order fixed; ppmSeq call order fixed) ---
  add_callset("Illumina")
  add_callset("Illumina_30x")
  add_callset("Illumina_100x")
  add_callset("Illumina_200x")

  add_callset("ppmSeq_multiread")
  add_callset("ppmSeq_putative_multiread")
  add_callset("ppmSeq_singleton_HC")

  add_callset("ppmSeq_multiread_30x")
  add_callset("ppmSeq_putative_multiread_30x")
  add_callset("ppmSeq_singleton_HC_30x")

  add_callset("ppmSeq_multiread_100x")
  add_callset("ppmSeq_putative_multiread_100x")
  add_callset("ppmSeq_singleton_HC_100x")

  add_callset("ppmSeq_multiread_200x")
  add_callset("ppmSeq_putative_multiread_200x")
  add_callset("ppmSeq_singleton_HC_200x")

  add_callset("UDSeq")
  add_callset("UDSeq_30x")

  add_callset("NanoSeq")

  add_callset("HiDEFseq_double_stranded")
  add_callset("HiDEFseq_single_stranded")

  # PTA callsets (01-16) come after PTA summaries (we will insert summaries right before PTA_01)
  for(i=1;i<=16;i++){
    p = sprintf("PTA_single_neuron_%02d", i)
    add_callset(p)
    pta_callset[i] = p
  }

  add_callset("AmpliSeq")

  # --- blacklist/annotation block (as you described 43-60) ---
  add_black("in_segdup")
  add_black("in_ucsc_gaps")
  add_black("in_satellites_rDNA")
  add_black("low_mappability")
  add_black("in_encode_blacklist")
  add_black("blacklist_any")
  add_black("blacklist_bitmask")
  add_black("blacklist_binary")
  add_black("gnomad_common")
  add_black("gnomad_af")
  add_black("in_known_indel_gatk")
  add_black("in_known_indel_mills1000G")
  add_black("in_known_indel_any")
  add_black("is_germline")
  add_black("in_repeatmasker")
  add_black("repeatmasker_class")
  add_black("repeatmasker_family")
  add_black("repeatmasker_name")

  # --- hg19 block (61-66) ---
  add_hg19("hg19_chr")
  add_hg19("hg19_pos")
  add_hg19("hg19_ref")
  add_hg19("hg19_alt")
  add_hg19("hg19_unmapped")
  add_hg19("hg19_pileup_unmapped")

  # --- germline-like flags (355-359) ---
  add_flag("is_germline_ILLUMINA_UPPER_CI")
  add_flag("is_germline_PPMSEQ_UPPER_CI")
  add_flag("is_germline_UDSEQ_UPPER_CI")
  add_flag("is_germline_NANOSEQ_UPPER_CI")
  add_flag("is_germline_any_bulk_upper_ci")
}

NR==1{
  for(i=1;i<=NF;i++){
    h = $i
    idx[h] = i
    header[i] = h
  }
  # required coords
  if(!idx["chr"] || !idx["pos"] || !idx["ref"] || !idx["alt"]) die("missing one of chr/pos/ref/alt")

  # verify callsets exist
  for(i=1;i<=arrN["CALLSET"];i++){
    nm = arr["CALLSET", i]
    if(!(nm in idx)) die("missing callset column: " nm)
    is_callset[idx[nm]] = 1
  }

  # verify black/hg19 exist
  for(i=1;i<=arrN["BLACK"];i++){
    nm = arr["BLACK", i]
    if(!(nm in idx)) die("missing blacklist/annot column: " nm)
    is_black[idx[nm]] = 1
  }
  for(i=1;i<=arrN["HG19"];i++){
    nm = arr["HG19", i]
    if(!(nm in idx)) die("missing hg19 column: " nm)
    is_hg19[idx[nm]] = 1
  }

  # flags: allow dropping final any-bulk
  for(i=1;i<=arrN["FLAGS"];i++){
    nm = arr["FLAGS", i]
    if(nm=="is_germline_any_bulk_upper_ci" && KEEP_ANY_BULK==0) continue
    if(!(nm in idx)) die("missing flag column: " nm)
    is_flag[idx[nm]] = 1
  }

  # collect BAM quant columns from header (67-354, but by pattern)
  for(i=1;i<=NF;i++){
    h = header[i]
    if(is_bam_quant(h)){
      bam_cols[++n_bam] = h
      is_bam[idx[h]] = 1
    }
  }
  if(n_bam==0) die("no BAM quant columns detected (pattern _REF_COUNT/_ALT_COUNT/... )")

  # sort BAM cols by category+method+pta+downsample (stable)
  for(i=1;i<=n_bam;i++){
    h = bam_cols[i]
    key = bam_sort_key(h)
    bam_key[h] = key
  }

  # build sorted unique list by key (gawk asort not guaranteed everywhere; implement via sorting keys lexically)
  for(i=1;i<=n_bam;i++){
    h = bam_cols[i]
    keys[++n_keys] = bam_key[h]
    key2col[bam_key[h]] = h
  }
  # simple lexical sort of keys (works because key encodes ranks)
  n_sorted = asort(keys, keys_sorted)

  # --- print header in final order ---
  # coords
  printf "%s%s%s%s%s%s%s", "chr", OFS, "pos", OFS, "ref", OFS, "alt"

  # callsets with call_ prefix + PTA summaries inserted before PTA_01
  inserted_pta_summary = 0
  for(i=1;i<=arrN["CALLSET"];i++){
    nm = arr["CALLSET", i]
    if(!inserted_pta_summary && nm=="PTA_single_neuron_01"){
      printf "%s%s%s%s%s", OFS, "call_PTA_number", OFS, "call_PTA_singleton", OFS
      printf "%s", "call_PTA_multiple"
      inserted_pta_summary = 1
    }
    printf "%s%s", OFS, "call_" nm
  }

  # blacklist/annotation
  for(i=1;i<=arrN["BLACK"];i++){
    nm = arr["BLACK", i]
    printf "%s%s", OFS, nm
  }

  # hg19
  for(i=1;i<=arrN["HG19"];i++){
    nm = arr["HG19", i]
    printf "%s%s", OFS, nm
  }

  # BAM quant (sorted)
  for(i=1;i<=n_sorted;i++){
    h = key2col[keys_sorted[i]]
    printf "%s%s", OFS, h
  }

  # flags
  for(i=1;i<=arrN["FLAGS"];i++){
    nm = arr["FLAGS", i]
    if(nm=="is_germline_any_bulk_upper_ci" && KEEP_ANY_BULK==0) continue
    printf "%s%s", OFS, nm
  }
  printf "\n"
  next
}

{
  # coords
  printf "%s%s%s%s%s%s%s", $idx["chr"], OFS, $idx["pos"], OFS, $idx["ref"], OFS, $idx["alt"]

  # PTA summary from PTA callset (NOT BAM quant)
  pta_n = 0
  for(i=1;i<=16;i++){
    v = $(idx[pta_callset[i]])
    if(v=="" ) v = 0
    if(v+0 != 0) pta_n++
  }
  pta_singleton = (pta_n==1 ? 1 : 0)
  pta_multiple  = (pta_n>=2 ? 1 : 0)

  # callsets with call_ prefix + insert PTA summaries
  inserted_pta_summary = 0
  for(i=1;i<=arrN["CALLSET"];i++){
    nm = arr["CALLSET", i]
    if(!inserted_pta_summary && nm=="PTA_single_neuron_01"){
      printf "%s%d%s%d%s%d", OFS, pta_n, OFS, pta_singleton, OFS, pta_multiple
      inserted_pta_summary = 1
    }
    printf "%s%s", OFS, $(idx[nm])
  }

  # blacklist/annotation
  for(i=1;i<=arrN["BLACK"];i++){
    nm = arr["BLACK", i]
    printf "%s%s", OFS, $(idx[nm])
  }

  # hg19
  for(i=1;i<=arrN["HG19"];i++){
    nm = arr["HG19", i]
    printf "%s%s", OFS, $(idx[nm])
  }

  # BAM quant (sorted)
  for(i=1;i<=n_sorted;i++){
    h = key2col[keys_sorted[i]]
    printf "%s%s", OFS, $(idx[h])
  }

  # flags
  for(i=1;i<=arrN["FLAGS"];i++){
    nm = arr["FLAGS", i]
    if(nm=="is_germline_any_bulk_upper_ci" && KEEP_ANY_BULK==0) continue
    printf "%s%s", OFS, $(idx[nm])
  }

  printf "\n"
}
' "$in_tsv" > "$out_tsv"

echo "[WROTE] $out_tsv" >&2
echo "[NOTE] any-bulk flag kept? $keep_any_bulk" >&2

