# Resolve portable paths from this configuration file, including under Slurm.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$CONFIG_DIR/../.." && pwd)"
IMAGE_PATH="${PACKAGE_ROOT}/containers/deepsomatic1.9.0_cpu.sif"
OUTPUT_DIR="${PACKAGE_ROOT}/results/deepsomatic"
REF_FASTA="${PACKAGE_ROOT}/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
TUMOR_BAM_PATH="${PACKAGE_ROOT}/inputs/7614_illumina_brain_somatic.bam"
NORMAL_BAM_PATH="${PACKAGE_ROOT}/inputs/7614_illumina_matchedNormal.bam"
BIND_PATHS="${PACKAGE_ROOT}"
