#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
mosdepth_blacklist_global.sh
  Create keep-regions = genome(bam contigs) minus blacklist, run mosdepth fast with flag filtering,
  and output plot-ready global coverage distributions (CDF/PDF) over the kept regions.

USAGE:
  mosdepth_blacklist_global.sh -b sample.bam -l blacklist.bed -o outdir -p prefix [options]

REQUIRED:
  -b  BAM (sorted + indexed .bai)
  -l  blacklist BED (0-based, half-open)
  -o  output directory
  -p  output prefix

OPTIONS:
  -t  threads for mosdepth decompression (default: 8)
  -q  min MAPQ (default: 20)
  -F  exclude-flag bits (default: 1796 = UNMAP+SECONDARY+QCFAIL+DUP)
  -T  comma thresholds for summary fractions (default: 1,5,10,20,30,50)
  -r  keep-contigs regex (default: '.*' i.e., all contigs in BAM)
  --keep-bed  path to write keep BED (default: outdir/prefix.keep.bed)
  --no-clean  keep temp files

OUTPUTS (in outdir):
  prefix.keep.bed                         # complement of blacklist within BAM contigs
  prefix.mosdepth.region.dist.txt          # cumulative dist over keep bed (total + per-contig)
  prefix.mosdepth.global.dist.txt          # whole-genome dist (regardless of keep bed)
  prefix.regions.bed.gz                    # mean depth per keep interval
  prefix.mosdepth.summary.txt              # mosdepth summary

  prefix.keep.global.cdf.tsv               # depth \t frac_ge (TOTAL over keep)
  prefix.keep.global.pdf.tsv               # depth \t frac_eq (TOTAL over keep)
  prefix.keep.global.summary.tsv           # key scalars + fractions at thresholds

DEPENDENCIES:
  mosdepth, samtools, bedtools, bgzip(optional), zcat/gzip
EOF
}

# defaults
THREADS=8
MAPQ=0
EXCL_FLAG=3840
THRESHOLDS="1,5,10,20,30,50"
CONTIG_RE=".*"
NO_CLEAN=0
KEEP_BED=""

# args
BAM=""
BLACKLIST=""
OUTDIR=""
PREFIX=""

# parse
while [[ $# -gt 0 ]]; do
  case "$1" in
    -b) BAM="$2"; shift 2;;
    -l) BLACKLIST="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    -p) PREFIX="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    -q) MAPQ="$2"; shift 2;;
    -F) EXCL_FLAG="$2"; shift 2;;
    -T) THRESHOLDS="$2"; shift 2;;
    -r) CONTIG_RE="$2"; shift 2;;
    --keep-bed) KEEP_BED="$2"; shift 2;;
    --no-clean) NO_CLEAN=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

if [[ -z "${BAM}" || -z "${BLACKLIST}" || -z "${OUTDIR}" || -z "${PREFIX}" ]]; then
  usage >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
if [[ -z "${KEEP_BED}" ]]; then
  KEEP_BED="${OUTDIR}/${PREFIX}.keep.bed"
fi

# sanity checks
command -v samtools >/dev/null
command -v bedtools >/dev/null
command -v mosdepth >/dev/null

[[ -f "${BAM}" ]] || { echo "Missing BAM: ${BAM}" >&2; exit 1; }
[[ -f "${BLACKLIST}" ]] || { echo "Missing blacklist BED: ${BLACKLIST}" >&2; exit 1; }
[[ -f "${BAM}.bai" || -f "${BAM%.bam}.bai" || -f "${BAM}.crai" || -f "${BAM%.cram}.crai" ]] || {
  echo "Missing alignment index (.bai or .crai). Please index: samtools index ${BAM}" >&2
  exit 1
}

tmpdir="$(mktemp -d "${OUTDIR}/${PREFIX}.tmp.XXXXXX")"
cleanup() { [[ "${NO_CLEAN}" -eq 1 ]] || rm -rf "${tmpdir}"; }
trap cleanup EXIT

GENOME="${tmpdir}/genome.sizes"
BL_FILT="${tmpdir}/blacklist.filtered.bed"
BL_SORT="${tmpdir}/blacklist.sorted.bed"

# 1) genome sizes from BAM index (contigs + lengths)
# idxstats columns: contig length mapped unmapped ; last row "*" is unmapped summary
samtools idxstats "${BAM}" \
  | awk -v re="${CONTIG_RE}" 'BEGIN{OFS="\t"} $1!="*" && $2>0 && $1 ~ re {print $1,$2}' \
  > "${GENOME}"

if [[ ! -s "${GENOME}" ]]; then
  echo "No contigs matched regex '${CONTIG_RE}' from BAM idxstats." >&2
  exit 1
fi

# 2) filter blacklist to contigs present in genome.sizes, then sort with that contig order
awk 'BEGIN{FS=OFS="\t"} NR==FNR{ok[$1]=1; next} ok[$1]' "${GENOME}" "${BLACKLIST}" \
  > "${BL_FILT}" || true

# If blacklist is empty after filtering, keep = whole genome
if [[ -s "${BL_FILT}" ]]; then
  bedtools sort -faidx "${GENOME}" -i "${BL_FILT}" > "${BL_SORT}"
  bedtools complement -i "${BL_SORT}" -g "${GENOME}" > "${KEEP_BED}"
else
  # whole genome as BED
  awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "${GENOME}" > "${KEEP_BED}"
fi

if [[ ! -s "${KEEP_BED}" ]]; then
  echo "Keep BED is empty (blacklist may cover everything?)." >&2
  exit 1
fi

# 3) mosdepth fast over keep-bed, excluding UNMAP+SECONDARY+QCFAIL+DUP by flags
OUTP="${OUTDIR}/${PREFIX}"
mosdepth -n --fast-mode -t "${THREADS}" -Q "${MAPQ}" -F "${EXCL_FLAG}" -f ./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa --by "${KEEP_BED}" \
  "${OUTP}" "${BAM}"


# 4) plot-ready global distributions over KEEP regions (use region.dist TOTAL line)
REGDIST="${OUTP}.mosdepth.region.dist.txt"
REGIONS="${OUTP}.regions.bed.gz"

CDF="${OUTDIR}/${PREFIX}.keep.global.cdf.tsv"
PDF="${OUTDIR}/${PREFIX}.keep.global.pdf.tsv"
SUM="${OUTDIR}/${PREFIX}.keep.global.summary.tsv"

# CDF: depth, frac_ge (from total rows)
awk 'BEGIN{OFS="\t"} $1=="total"{print $2,$3}' "${REGDIST}" > "${CDF}"

# PDF: depth, frac_eq = frac_ge(d) - frac_ge(d+1)
awk 'BEGIN{OFS="\t"}
     $1=="total"{
       d=$2; p=$3;
       if(seen){
         print prev_d, (prev_p - p);
       }
       prev_d=d; prev_p=p; seen=1;
     }
     END{
       if(seen){ print prev_d, prev_p; }
     }' "${REGDIST}" > "${PDF}"

# Weighted mean depth over KEEP regions from regions.bed.gz
# regions.bed.gz columns: chr start end mean (because we wrote 3-col keep bed)
MEAN_DEPTH=$(zcat "${REGIONS}" | awk 'BEGIN{FS=OFS="\t"} {len=$3-$2; sum+=len*$4; L+=len} END{if(L>0) printf("%.6f", sum/L); else print "nan"}')
KEEP_BP=$(awk 'BEGIN{FS=OFS="\t"} {L+=($3-$2)} END{print L+0}' "${KEEP_BED}")

# Fractions at thresholds from CDF (if threshold > max depth, return 0)
# We also capture depth=0 fraction>=0 (should be 1.0)
{
  echo -e "metric\tvalue"
  echo -e "bam\t${BAM}"
  echo -e "blacklist_bed\t${BLACKLIST}"
  echo -e "keep_bed\t${KEEP_BED}"
  echo -e "exclude_flag\t${EXCL_FLAG}"
  echo -e "min_mapq\t${MAPQ}"
  echo -e "fast_mode\t1"
  echo -e "keep_bases\t${KEEP_BP}"
  echo -e "mean_depth_keep\t${MEAN_DEPTH}"
  echo -e "dist_file\t${REGDIST}"
  echo -e "cdf_file\t${CDF}"
  echo -e "pdf_file\t${PDF}"
} > "${SUM}"

# add threshold rows
IFS=',' read -r -a THR <<< "${THRESHOLDS}"
for t in "${THR[@]}"; do
  frac=$(awk -v T="${t}" 'BEGIN{FS=OFS="\t"; f=0}
                         $1=="total" && $2==T {f=$3}
                         END{printf("%.6f", f)}' "${REGDIST}")
  echo -e "frac_ge_${t}x\t${frac}" >> "${SUM}"
done

echo "Done."
echo "Keep regions: ${KEEP_BED}"
echo "Plot-ready:   ${CDF} (CDF), ${PDF} (PDF)"
echo "Summary:      ${SUM}"