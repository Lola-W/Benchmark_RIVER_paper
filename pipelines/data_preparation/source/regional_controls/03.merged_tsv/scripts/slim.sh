#!/usr/bin/env bash
set -euo pipefail

in="benchmark_with_downsampling_quant_added.flags.tsv"
out="slim.benchmark_with_downsampling_quant_added.flags.with_pta_which.tsv"

awk -F'\t' '
BEGIN{
  OFS="\t"
  uci_thr=0.52
  alt_min=1
  depth_min=3
}

NR==1{
  # Cache header names and detect indices
  for(i=1;i<=NF;i++){
    h[i]=$i

    # Indices for PTA per-neuron pileup metrics (used for computing which-lists)
    if(match(h[i], /^PTA_SINGLE_NEURON_([0-9][0-9])_ALT_COUNT$/, m)){
      n = m[1] + 0
      pta_alt_idx[n]=i
    }
    if(match(h[i], /^PTA_SINGLE_NEURON_([0-9][0-9])_REF_COUNT$/, m)){
      n = m[1] + 0
      pta_ref_idx[n]=i
    }
    if(match(h[i], /^PTA_SINGLE_NEURON_([0-9][0-9])_UPPER_CI$/, m)){
      n = m[1] + 0
      pta_uci_idx[n]=i
    }

    # Where to insert new columns
    if(h[i] == "is_maf_upper_ci_ge_052_PTA_number") maf_pta_number_i=i
    if(h[i] == "is_alt_count_ge_1_PTA_number")     alt_pta_number_i=i
  }

  # Print slim header with insertions, while dropping unwanted columns
  outc=0
  for(i=1;i<=NF;i++){
    col=h[i]

    drop=0

    # Drop VCF-level PTA individual calls
    if(col ~ /^call_PTA_single_neuron_/) drop=1

    # Drop pileup-level PTA single neuron metrics (we will compute which-list before dropping)
    if(col ~ /^PTA_SINGLE_NEURON_[0-9][0-9]_(ALT_COUNT|IS_GREATER|LOWER_CI|MAF|REF_COUNT|UPPER_CI)$/) drop=1

    # Drop unrelated / negative control pileups (all metrics)
    if(col ~ /^HG19_ILLUMINA_7669_/) drop=1
    if(col ~ /^ILLUMINA_6566_BLOOD_/) drop=1
    if(col ~ /^ILLUMINA_6566_SPERM_/) drop=1
    if(col ~ /^PPMSEQ_6566_/) drop=1
    if(col ~ /^UDSEQ_6566_/) drop=1
    if(col ~ /^NANOSEQ_6566_/) drop=1
    if(col ~ /^HG19_PTA_6566_1_/) drop=1
    if(col ~ /^HG19_PTA_6566_2_/) drop=1
    if(col ~ /^HG19_PTA_8420_1_/) drop=1
    if(col ~ /^HG19_PTA_8420_2_/) drop=1

    if(!drop){
      # Insert new header field right before PTA_number summaries
      if(i==maf_pta_number_i){
        out[++outc]="is_maf_upper_ci_ge_052_PTA_which"
      }
      if(i==alt_pta_number_i){
        out[++outc]="is_alt_count_ge_1_PTA_which"
      }
      out[++outc]=col
    }
  }

  for(k=1;k<=outc;k++){
    printf "%s%s", out[k], (k==outc?ORS:OFS)
  }
  next
}

{
  # Compute which-lists from per-neuron pileups (even though they will be dropped)
  maf_list=""
  alt_list=""

  for(n=1;n<=16;n++){
    ai=pta_alt_idx[n]; ri=pta_ref_idx[n]; ui=pta_uci_idx[n]

    alt=(ai>0 && $(ai)!="" ? $(ai)+0 : 0)
    ref=(ri>0 && $(ri)!="" ? $(ri)+0 : 0)
    uci_str=(ui>0 ? $(ui) : "")

    # alt>=1 list
    if(alt >= alt_min){
      if(alt_list=="") alt_list=n
      else alt_list=alt_list "," n
    }

    # maf_upper_ci list:
    # ALT>=1 AND ALT+REF>=3 AND UPPER_CI not NA AND UPPER_CI>0.52
    if(alt >= alt_min && (alt+ref) >= depth_min){
      if(uci_str != "" && uci_str != "NA" && uci_str != "."){
        uci=uci_str+0
        if(uci > uci_thr){
          if(maf_list=="") maf_list=n
          else maf_list=maf_list "," n
        }
      }
    }
  }

  # Build slim row with insertions + drops
  outc=0
  for(i=1;i<=NF;i++){
    col=h[i]
    drop=0

    if(col ~ /^call_PTA_single_neuron_/) drop=1
    if(col ~ /^PTA_SINGLE_NEURON_[0-9][0-9]_(ALT_COUNT|IS_GREATER|LOWER_CI|MAF|REF_COUNT|UPPER_CI)$/) drop=1

    if(col ~ /^HG19_ILLUMINA_7669_/) drop=1
    if(col ~ /^ILLUMINA_6566_BLOOD_/) drop=1
    if(col ~ /^ILLUMINA_6566_SPERM_/) drop=1
    if(col ~ /^PPMSEQ_6566_/) drop=1
    if(col ~ /^UDSEQ_6566_/) drop=1
    if(col ~ /^NANOSEQ_6566_/) drop=1
    if(col ~ /^HG19_PTA_6566_1_/) drop=1
    if(col ~ /^HG19_PTA_6566_2_/) drop=1
    if(col ~ /^HG19_PTA_8420_1_/) drop=1
    if(col ~ /^HG19_PTA_8420_2_/) drop=1

    if(!drop){
      if(i==maf_pta_number_i){
        out[++outc]=maf_list
      }
      if(i==alt_pta_number_i){
        out[++outc]=alt_list
      }
      out[++outc]=$i
    }
  }

  for(k=1;k<=outc;k++){
    printf "%s%s", out[k], (k==outc?ORS:OFS)
  }
}
' "$in" > "$out"

echo "[DONE] $out"
ls -lh "$in" "$out" | sed "s/^/[SIZE] /"

