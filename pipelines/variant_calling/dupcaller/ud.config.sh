# Resolve portable paths from this configuration file, including under Slurm.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$CONFIG_DIR/.." && pwd)"
REF="${PACKAGE_ROOT}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
NOISE="${PACKAGE_ROOT}/resources/hg38.noise.bed.gz"
GERMLINE="${PACKAGE_ROOT}/inputs/41411195_F_illumina_matchedNormal.haplotypecaller.SNV_INDEL.VQSR_INFERRED.vcf.gz"
TUMOR="${PACKAGE_ROOT}/inputs/7614-Cortex-L-T_D.bam"
NORMAL="${PACKAGE_ROOT}/inputs/7614-Kidney-L_D.bam"
OUTPREFIX="${PACKAGE_ROOT}/results/ud/7614-Cortex-L-T_D"
