#!/usr/bin/env bash
#SBATCH --job-name=trinuc_round2
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=7-00:00:00
#SBATCH --array=0-22
#SBATCH --output=log/%x.%A_%a.out
#SBATCH --error=log/%x.%A_%a.err
#SBATCH --export=ALL

set -euo pipefail
set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

THREADS=12
MIN_BQ=13
EXCL_FLAGS=3844

CHR_LIST=()
for c in {1..22}; do CHR_LIST+=("chr$c"); done
CHR_LIST+=("chrX")
CHR="${CHR_LIST[$SLURM_ARRAY_TASK_ID]}"

TMPROOT="${OUTROOT}/round2_tmp/${KEY}"
mkdir -p "${TMPROOT}"

LEN=$(awk -v C="${CHR}" '$1==C{print $2; exit}' "${FA}.fai")

echo -e "${CHR}\t0\t${LEN}" > "${TMPROOT}/${CHR}.full.bed"
awk -v C="${CHR}" '$1==C' "${BLACK}" \
| sort -k1,1 -k2,2n -k3,3n > "${TMPROOT}/${CHR}.black.bed"

bedtools subtract -sorted \
  -a "${TMPROOT}/${CHR}.full.bed" \
  -b "${TMPROOT}/${CHR}.black.bed" \
  > "${TMPROOT}/${CHR}.eff.bed"

samtools depth -@ ${THREADS} -q ${MIN_BQ} \
  -b "${TMPROOT}/${CHR}.eff.bed" \
  -r "${CHR}" \
  -G ${EXCL_FLAGS} \
  "${BAM}" \
| awk -v L=${LEN} '{pos=$2; if(pos<2 || pos>L-1) next; print pos"\t"$3}' \
> "${TMPROOT}/${CHR}.pos_depth.tsv"

awk -v C=${CHR} '{s=$1-1; e=$1+1; print C":"s"-"e}' \
"${TMPROOT}/${CHR}.pos_depth.tsv" > "${TMPROOT}/${CHR}.regions.txt"

samtools faidx "${FA}" -r "${TMPROOT}/${CHR}.regions.txt" \
| awk '
function rc3(s,   i,b,r){
r=""; for(i=3;i>=1;i--){
b=substr(s,i,1);
if(b=="A") r=r "T";
else if(b=="T") r=r "A";
else if(b=="C") r=r "G";
else if(b=="G") r=r "C";
else r=r "N";
}
return r;
}
/^>/{
gsub(/^>/,""); split($0,a,":"); split(a[2],b,"-"); pos=b[1]+1;
getline seq; gsub(/[ \r\n\t]/,"",seq);
seq=toupper(seq);
if(length(seq)!=3) seq="NNN";
c=substr(seq,2,1);
if(c=="A"||c=="G") seq=rc3(seq);
print pos"\t"seq
}' > "${TMPROOT}/${CHR}.tri.tsv"

paste "${TMPROOT}/${CHR}.pos_depth.tsv" "${TMPROOT}/${CHR}.tri.tsv" \
| awk -v C=${CHR} '{start=$1-1; print C"\t"start"\t"$1"\t"$2"\t"$4}' \
> "${TMPROOT}/${CHR}.tsv"
