#!/bin/bash
# Stage a run directory so `nucleosome-adaptive-sampling` can find what it
# reads as relative paths: ./ref/<linker_length>/<linker_length>rep.pdb,
# ./mdp_features.dat, ./traj_dat/*.xtc
#
# Usage:
#   scripts/setup_adaptive_run.sh <run_dir> <linker_length>
#
# Example:
#   scripts/setup_adaptive_run.sh /scratch/me/run1 15
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="${1:?usage: setup_adaptive_run.sh <run_dir> <linker_length>}"
LINKER_LENGTH="${2:?usage: setup_adaptive_run.sh <run_dir> <linker_length>}"

SRC_REF_DIR="${REPO_ROOT}/ref/${LINKER_LENGTH}"
if [[ ! -d "${SRC_REF_DIR}" ]]; then
    echo "ERROR: no reference structures for linker length ${LINKER_LENGTH} (expected ${SRC_REF_DIR})" >&2
    exit 1
fi

DST_REF_DIR="${RUN_DIR}/ref/${LINKER_LENGTH}"
mkdir -p "${DST_REF_DIR}" "${RUN_DIR}/traj_dat" "${RUN_DIR}/CV"

# adaptive.py derives topology_file as ref/<linker_length>/<linker_length>rep.pdb
# by default (override with --topology-file if your layout differs).
ln -sf "${SRC_REF_DIR}/${LINKER_LENGTH}rep.pdb" "${DST_REF_DIR}/${LINKER_LENGTH}rep.pdb"
ln -sf "${SRC_REF_DIR}/DNA.pdb" "${DST_REF_DIR}/DNA.pdb"
ln -sf "${SRC_REF_DIR}/core_all.pdb" "${DST_REF_DIR}/core_all.pdb"

if [[ -z "$(ls -A "${RUN_DIR}/traj_dat" 2>/dev/null)" ]]; then
    echo "NOTE: ${RUN_DIR}/traj_dat is empty. Put the original seed trajectories" >&2
    echo "      (*.xtc, n_original_sims=320 by default) there before running." >&2
fi

nucleosome-feature-writer \
    --pdb "${DST_REF_DIR}/${LINKER_LENGTH}rep.pdb" \
    --linker-length "${LINKER_LENGTH}" \
    --output "${RUN_DIR}/mdp_features.dat"

python3 - <<'PYEOF'
import openmm as mm
names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
if "CUDA" not in names:
    print(f"WARNING: no CUDA platform available (found: {names}). "
          "run_one_adaptive_simulation() will fail for every seed.")
PYEOF

echo "Staged ${RUN_DIR}:"
echo "  ${DST_REF_DIR}/${LINKER_LENGTH}rep.pdb -> ${SRC_REF_DIR}/${LINKER_LENGTH}rep.pdb"
echo "  ${RUN_DIR}/mdp_features.dat"
echo "  ${RUN_DIR}/traj_dat/  (populate with *.xtc before running)"
echo "  ${RUN_DIR}/CV/        (will be filled in by the pipeline)"
echo
echo "Next: cd ${RUN_DIR} && nucleosome-adaptive-sampling --n-rounds 8 --n-seed-per-round 16 --linker-length ${LINKER_LENGTH}"
