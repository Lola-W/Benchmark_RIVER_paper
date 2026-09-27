#!/usr/bin/env bash
set -euo pipefail

MASTER="benchmark_snv_matrix.annot.with_hg19.plus_pileup.tsv"
BASE_DIR="down_baseQ13_hg38/result"
OUT="benchmark_with_downsampling_quant_added.tsv"

TMPDIR=$(mktemp -d)

################################
# 1. Sort MASTER
################################
echo "[1] Sort MASTER"
(head -n1 "$MASTER" && tail -n+2 "$MASTER" \
    | sort -k1,1V -k2,2n -k3,3 -k4,4) > "$TMPDIR/master.sorted.tsv"

################################
# 2. Prepare downsample blocks
################################
prep_ds () {
    FILE="$1"
    LABEL="$2"

    TMP="$TMPDIR/${LABEL}.sorted.tsv"

    (head -n1 "$FILE" && tail -n+2 "$FILE" \
        | sort -k1,1V -k2,2n -k3,3 -k4,4) > "$TMP"

    # Coordinate check
    diff \
      <(awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' "$TMPDIR/master.sorted.tsv") \
      <(awk 'NR>1{print $1"\t"$2"\t"$3"\t"$4}' "$TMP") \
      > /dev/null || { echo "Coordinate mismatch in $LABEL"; exit 1; }

    # Block 생성 (6 columns)
    awk -F'\t' -v L="$LABEL" -v OFS='\t' '
    NR==1{
        print L"_REF_COUNT",L"_ALT_COUNT",L"_MAF",L"_LOWER_CI",L"_UPPER_CI",L"_IS_GREATER"
        next
    }
    { print $5,$6,$7,$8,$9,$10 }
    ' "$TMP" > "$TMPDIR/${LABEL}.block"
}

echo "[2] Preparing blocks"

prep_ds "$BASE_DIR/illumina_30x/baseQ13.big_with_downsampling.maf_ci.tsv"   ILLUMINA_30x
prep_ds "$BASE_DIR/illumina_100x/baseQ13.big_with_downsampling.maf_ci.tsv"  ILLUMINA_100x
prep_ds "$BASE_DIR/illumina_200x/baseQ13.big_with_downsampling.maf_ci.tsv"  ILLUMINA_200x

prep_ds "$BASE_DIR/ppmseq_30x/baseQ13.big_with_downsampling.maf_ci.tsv"     PPMSEQ_30x
prep_ds "$BASE_DIR/ppmseq_100x/baseQ13.big_with_downsampling.maf_ci.tsv"    PPMSEQ_100x
prep_ds "$BASE_DIR/ppmseq_200x/baseQ13.big_with_downsampling.maf_ci.tsv"    PPMSEQ_200x

prep_ds "$BASE_DIR/udseq_30x/baseQ13.big_with_downsampling.maf_ci.tsv"      UDSEQ_30x

################################
# 3. Structural insertion
################################
echo "[3] Structural insertion"

paste \
    "$TMPDIR/master.sorted.tsv" \
    "$TMPDIR/ILLUMINA_30x.block" \
    "$TMPDIR/ILLUMINA_100x.block" \
    "$TMPDIR/ILLUMINA_200x.block" \
    "$TMPDIR/PPMSEQ_30x.block" \
    "$TMPDIR/PPMSEQ_100x.block" \
    "$TMPDIR/PPMSEQ_200x.block" \
    "$TMPDIR/UDSEQ_30x.block" \
> "$TMPDIR/full_paste.tmp"

awk -F'\t' -v OFS='\t' '
NR==1{
    # find anchor positions in original master
    for(i=1;i<=NF;i++){
        if($i=="ILLUMINA_IS_GREATER") illum_end=i
        if($i=="PPMSEQ_IS_GREATER") ppm_end=i
        if($i=="UDSEQ_IS_GREATER") ud_end=i
    }

    orig_nf=NF-42   # original column count
}

{
    # split full line
    split($0,a,"\t")

    # original columns
    o_nf=orig_nf

    # appended block offset
    offset=o_nf

    # illuminate insertion indices
    illum_ins=illum_end
    ppm_ins=ppm_end
    ud_ins=ud_end

    # build new line
    out=""

    for(i=1;i<=o_nf;i++){

        out = (i==1 ? a[i] : out OFS a[i])

        if(i==illum_ins){
            for(j=1;j<=18;j++)
                out=out OFS a[o_nf+j]
        }

        if(i==ppm_ins){
            for(j=19;j<=36;j++)
                out=out OFS a[o_nf+j]
        }

        if(i==ud_ins){
            for(j=37;j<=42;j++)
                out=out OFS a[o_nf+j]
        }
    }

    print out
}
' "$TMPDIR/full_paste.tmp" > "$OUT"

################################
# 4. Final check
################################
echo "[4] Row count"
wc -l "$OUT"

rm -rf "$TMPDIR"

echo "DONE"

