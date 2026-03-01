#!/usr/bin/env bash
set -euo pipefail

# Step 2: XYZ output -> Br/Btheta/Bphi output
#
# Usage:
#   ./run_xyz_to_brtp.sh <input_xyz_file> [output_brtp_file]
#
# This step is deterministic post-processing:
# - no source physics, only coordinate-basis rotation

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

XYZ_IN="${1:-}"
BRTP_OUT="${2:-}"

if [[ -z "${XYZ_IN}" ]]; then
  echo "Usage: $0 <input_xyz_file> [output_brtp_file]"
  exit 2
fi

if [[ ! -f "${XYZ_IN}" ]]; then
  echo "ERROR: xyz input file not found: ${XYZ_IN}"
  exit 2
fi

if [[ -z "${BRTP_OUT}" ]]; then
  base="$(basename "${XYZ_IN}")"
  base="${base%_xyz.txt}"
  if [[ "${base}" == "$(basename "${XYZ_IN}")" ]]; then
    base="${base%.txt}"
  fi
  mkdir -p "${ROOT_DIR}/output"
  BRTP_OUT="${ROOT_DIR}/output/${base}_brtp.txt"
fi

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  # converter is tiny, but we keep build behavior consistent with other wrappers.
  "${ROOT_DIR}/build_gravmag_tools.sh"
fi

echo "+ step2 convert: ./gravmag_xyz_to_brtp ${XYZ_IN} ${BRTP_OUT}"
"${ROOT_DIR}/gravmag_xyz_to_brtp" "${XYZ_IN}" "${BRTP_OUT}"

echo "Done: ${BRTP_OUT}"
