# =========================
# File: submit_tri32_covhist_array.sh
# =========================
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="./work/regional_controls/06.trinuc_composition"
cd "${BASE_DIR}"

mkdir -p log

sbatch trinuc32_and_covhist_array.sbatch
