#!/usr/bin/env bash
# Build the Lace 2.0 Apptainer container.
# Run from the Lace repo root:
#   cd /path/to/Lace && bash build.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEF="${SCRIPT_DIR}/Lace_v2.def"
SIF="${SCRIPT_DIR}/lace_2.0.sif"

echo "Building ${SIF} from ${DEF} ..."
apptainer build --fakeroot "${SIF}" "${DEF}"
echo "Done. Test with:  apptainer exec ${SIF} Lace --help"
