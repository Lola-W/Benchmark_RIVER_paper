#!/usr/bin/env bash
#SBATCH --job-name=duplex_cat_nano_6566
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
# Set site-specific scheduler options in the submitting environment.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --export=ALL
#SBATCH --output=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/5hidef_duplex_cat.%j.log
#SBATCH --error=./work/regional_controls/05.cleaned_bam/tmp/_CHR_LOGS/5hidef_duplex_cat.%j.err

set -euo pipefail

set +u
# Activate the required software environment before submission.
set -u
# Activate the required software environment before submission.

BASE=./work/regional_controls/05.cleaned_bam
cd "$BASE"

OUTDIR="$BASE"
mkdir -p "$OUTDIR"

OUT="$OUTDIR/duplex_clean.NanoSeq_6566.chr1-22_X.bam"
TMPDIR="$BASE/tmp/_cat_tmp"
mkdir -p "$TMPDIR"

# chr1..chr22, chrX (V-sort order)
BAMS=()
for c in {1..22} X; do
  f="$BASE/duplex_clean.S6566_Hpy.chr${c}.bam"
  [[ -s "$f" ]] || { echo "Missing: $f" >&2; exit 1; }
  BAMS+=("$f")
done

# Quick header @SQ compatibility check
REFHDR="$TMPDIR/ref.sq"
samtools view -H "${BAMS[0]}" | grep '^@SQ' > "$REFHDR"
for f in "${BAMS[@]:1}"; do
  samtools view -H "$f" | grep '^@SQ' | diff -q "$REFHDR" - >/dev/null || {
    echo "Header @SQ mismatch: $f" >&2
    exit 1
  }
done

# cat: single-thread (no -@ anyway)
samtools cat -h <(samtools view -H "${BAMS[0]}") -o "$OUT" "${BAMS[@]}"

# index: 12 threads
samtools index -@ "${SLURM_CPUS_PER_TASK}" "$OUT"

samtools quickcheck -v "$OUT"
samtools idxstats "$OUT" | head
echo "DONE: $OUT"
