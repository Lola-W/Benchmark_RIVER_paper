#!/usr/bin/env bash
set -euo pipefail

tsv_in="slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.tsv"

vcf_mr_30="ppm_vcf_after_qc/step07-1_trinuc_denoised_multiread_30x.nonMISMATCH.vcf"
vcf_mr_100="ppm_vcf_after_qc/step07-1_trinuc_denoised_multiread_100x.nonMISMATCH.vcf"
vcf_mr_200="ppm_vcf_after_qc/step07-1_trinuc_denoised_multiread_200x.nonMISMATCH.vcf"
vcf_mr_300="ppm_vcf_after_qc/step07-1_trinuc_denoised_multiread_300x.nonMISMATCH.vcf"

vcf_pmr_30="ppm_vcf_after_qc/step09-2_trinuc_denoised_singleton_putative_multiread_30x.dupraw_gt1_eq_rawmulti.vcf"
vcf_pmr_100="ppm_vcf_after_qc/step09-2_trinuc_denoised_singleton_putative_multiread_100x.dupraw_gt1_eq_rawmulti.vcf"
vcf_pmr_200="ppm_vcf_after_qc/step09-2_trinuc_denoised_singleton_putative_multiread_200x.dupraw_gt1_eq_rawmulti.vcf"
vcf_pmr_300="ppm_vcf_after_qc/step09-2_trinuc_denoised_singleton_putative_multiread_300x.dupraw_gt1_eq_rawmulti.vcf"

tsv_out="postqc.slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.tsv"
tsv_out_drop="postqc.slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.drop_allzero_call.tsv"

for f in "$tsv_in" \
  "$vcf_mr_30" "$vcf_mr_100" "$vcf_mr_200" "$vcf_mr_300" \
  "$vcf_pmr_30" "$vcf_pmr_100" "$vcf_pmr_200" "$vcf_pmr_300"
do
  [[ -f "$f" ]] || { echo "[ERROR] missing: $f" >&2; exit 2; }
done

gawk -F'\t' -v OFS='\t' \
  -v OUT1="$tsv_out" -v OUT2="$tsv_out_drop" \
  -v MR30="$vcf_mr_30" -v MR100="$vcf_mr_100" -v MR200="$vcf_mr_200" -v MR300="$vcf_mr_300" \
  -v PMR30="$vcf_pmr_30" -v PMR100="$vcf_pmr_100" -v PMR200="$vcf_pmr_200" -v PMR300="$vcf_pmr_300" \
'
# -------------------------
# Build QC membership sets from VCFs (key = chr\tpos\tref\talt)
# -------------------------
function add_vcf_key(map, chr, pos, ref, alt,    k){
  k = chr "\t" pos "\t" ref "\t" alt
  map[k] = 1
}

# Read VCFs first (ARGV order is controlled by the command line below)
FNR==1{
  # no-op: keep
}

# Identify which file is being read (gawk ARGIND)
ARGIND==1 { # MR30
  if($0 ~ /^#/) next
  add_vcf_key(mr30, $1, $2, $4, $5)
  next
}
ARGIND==2 { # MR100
  if($0 ~ /^#/) next
  add_vcf_key(mr100, $1, $2, $4, $5)
  next
}
ARGIND==3 { # MR200
  if($0 ~ /^#/) next
  add_vcf_key(mr200, $1, $2, $4, $5)
  next
}
ARGIND==4 { # MR300
  if($0 ~ /^#/) next
  add_vcf_key(mr300, $1, $2, $4, $5)
  next
}
ARGIND==5 { # PMR30
  if($0 ~ /^#/) next
  add_vcf_key(pmr30, $1, $2, $4, $5)
  next
}
ARGIND==6 { # PMR100
  if($0 ~ /^#/) next
  add_vcf_key(pmr100, $1, $2, $4, $5)
  next
}
ARGIND==7 { # PMR200
  if($0 ~ /^#/) next
  add_vcf_key(pmr200, $1, $2, $4, $5)
  next
}
ARGIND==8 { # PMR300
  if($0 ~ /^#/) next
  add_vcf_key(pmr300, $1, $2, $4, $5)
  next
}

# -------------------------
# TSV starts at ARGIND==9
# - Overwrite ppmSeq call columns using QC membership sets
# - Output:
#   OUT1: same rows, updated calls
#   OUT2: drop rows where all call_* are 0 (after overwrite)
# - Print before/after counts to stderr
# -------------------------
ARGIND==9 && NR==FNR{
  # never reached; kept for clarity
}

ARGIND==9 && FNR==1{
  # cache header + locate indices
  for(i=1;i<=NF;i++){
    h[i]=$i

    # collect all call_* columns for all-zero filtering
    if(h[i] ~ /^call_/) call_idx[++n_call]=i

    # ppmSeq call columns (NOTE: 300x has NO suffix in TSV)
    if(h[i]=="call_ppmSeq_multiread_30x") mr30_i=i
    if(h[i]=="call_ppmSeq_multiread_100x") mr100_i=i
    if(h[i]=="call_ppmSeq_multiread_200x") mr200_i=i
    if(h[i]=="call_ppmSeq_multiread") mr300_i=i

    if(h[i]=="call_ppmSeq_putative_multiread_30x") pmr30_i=i
    if(h[i]=="call_ppmSeq_putative_multiread_100x") pmr100_i=i
    if(h[i]=="call_ppmSeq_putative_multiread_200x") pmr200_i=i
    if(h[i]=="call_ppmSeq_putative_multiread") pmr300_i=i

    if(h[i]=="chr") chr_i=i
    if(h[i]=="pos") pos_i=i
    if(h[i]=="ref") ref_i=i
    if(h[i]=="alt") alt_i=i
  }

  # minimal sanity
  if(!(mr30_i && mr100_i && mr200_i && mr300_i && pmr30_i && pmr100_i && pmr200_i && pmr300_i)){
    print "[ERROR] missing expected ppmSeq call columns in TSV header" > "/dev/stderr"
    exit 3
  }

  # write header to both outputs
  print $0 > OUT1
  print $0 > OUT2
  next
}

ARGIND==9 && FNR>1{
  # count BEFORE overwrite (from existing TSV)
  if($(mr30_i)+0==1) b_mr30++
  if($(mr100_i)+0==1) b_mr100++
  if($(mr200_i)+0==1) b_mr200++
  if($(mr300_i)+0==1) b_mr300++
  if($(pmr30_i)+0==1) b_pmr30++
  if($(pmr100_i)+0==1) b_pmr100++
  if($(pmr200_i)+0==1) b_pmr200++
  if($(pmr300_i)+0==1) b_pmr300++

  # build key from TSV coordinates
  k = $(chr_i) "\t" $(pos_i) "\t" $(ref_i) "\t" $(alt_i)

  # overwrite ppmSeq calls using QC VCF membership sets
  $(mr30_i)  = (k in mr30  ? 1 : 0)
  $(mr100_i) = (k in mr100 ? 1 : 0)
  $(mr200_i) = (k in mr200 ? 1 : 0)
  $(mr300_i) = (k in mr300 ? 1 : 0)

  $(pmr30_i)  = (k in pmr30  ? 1 : 0)
  $(pmr100_i) = (k in pmr100 ? 1 : 0)
  $(pmr200_i) = (k in pmr200 ? 1 : 0)
  $(pmr300_i) = (k in pmr300 ? 1 : 0)

  # count AFTER overwrite
  if($(mr30_i)+0==1) a_mr30++
  if($(mr100_i)+0==1) a_mr100++
  if($(mr200_i)+0==1) a_mr200++
  if($(mr300_i)+0==1) a_mr300++
  if($(pmr30_i)+0==1) a_pmr30++
  if($(pmr100_i)+0==1) a_pmr100++
  if($(pmr200_i)+0==1) a_pmr200++
  if($(pmr300_i)+0==1) a_pmr300++

  # write full (updated) row
  print $0 > OUT1

  # drop rows where ALL call_* == 0 (after overwrite)
  any=0
  for(j=1;j<=n_call;j++){
    idx = call_idx[j]
    if($(idx)+0==1){ any=1; break }
  }
  if(any){
    print $0 > OUT2
  } else {
    dropped++
  }
  total++
  next
}

END{
  # report summary to stderr
  print "==== ppmSeq call counts (TSV) BEFORE -> AFTER QC overwrite ====" > "/dev/stderr"
  print "MR  30x", b_mr30+0, "->", a_mr30+0 > "/dev/stderr"
  print "MR 100x", b_mr100+0, "->", a_mr100+0 > "/dev/stderr"
  print "MR 200x", b_mr200+0, "->", a_mr200+0 > "/dev/stderr"
  print "MR 300x", b_mr300+0, "->", a_mr300+0 > "/dev/stderr"
  print "PMR  30x", b_pmr30+0, "->", a_pmr30+0 > "/dev/stderr"
  print "PMR 100x", b_pmr100+0, "->", a_pmr100+0 > "/dev/stderr"
  print "PMR 200x", b_pmr200+0, "->", a_pmr200+0 > "/dev/stderr"
  print "PMR 300x", b_pmr300+0, "->", a_pmr300+0 > "/dev/stderr"
  print "Rows total:", total+0, " rows dropped(all call_*==0):", dropped+0 > "/dev/stderr"

  print "[DONE] wrote:", OUT1 > "/dev/stderr"
  print "[DONE] wrote:", OUT2 > "/dev/stderr"
}
' \
  "$vcf_mr_30" "$vcf_mr_100" "$vcf_mr_200" "$vcf_mr_300" \
  "$vcf_pmr_30" "$vcf_pmr_100" "$vcf_pmr_200" "$vcf_pmr_300" \
  "$tsv_in"

echo "[OUT] $tsv_out"
echo "[OUT] $tsv_out_drop"

