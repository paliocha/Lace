#!/usr/bin/env bash
# Build the lace-fast 2.0 Singularity/Apptainer container.
#
# Run from the Lace repo root:
#   cd /path/to/Lace && bash lace-fast/build.sh
#
# The lace-fast/ source tree is copied into the container at build time
# via the %files section in Lace_fast.def.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="${SCRIPT_DIR}/Lace_fast.def"
SIF_FILE="${SCRIPT_DIR}/lace_fast_2.0.sif"

echo "Building lace-fast container..."
echo "  Definition: ${DEF_FILE}"
echo "  Output:     ${SIF_FILE}"

# Build from repo root so %files paths resolve correctly
cd "${SCRIPT_DIR}/.."

apptainer build --fakeroot "${SIF_FILE}" "${DEF_FILE}"

echo "Done: ${SIF_FILE}"
echo "Test with: apptainer exec ${SIF_FILE} lace-fast --help"
