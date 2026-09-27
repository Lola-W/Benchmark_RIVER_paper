# Resolve portable paths from this configuration file, including under Slurm.
CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(cd "$CONFIG_DIR/../.." && pwd)"
IMAGE_PATH="${PACKAGE_ROOT}/containers/deepmosaic_latest.sif"
OUTPUT_DIR="${PACKAGE_ROOT}/results/deepmosaic"
INPUT_TSV="${PACKAGE_ROOT}/illumina/deepmosaic/samples.tsv"
ANNOVAR_DIR="${PACKAGE_ROOT}/software/annovar"
BIND_PATHS="${PACKAGE_ROOT}"
