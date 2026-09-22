#!/usr/bin/env bash
set -euo pipefail

# XYZ output -> Br/Btheta/Bphi output
#
# Usage:
#   ./run_xyz_to_brtp.sh <input_xyz_file> [output_brtp_file]
#
# this is strictly post-processing, only coordinate-basis rotation

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
  # keep derived inputs in their owning case's output directory
  output_dir="$(cd -- "$(dirname -- "${XYZ_IN}")" && pwd)"
  parent_dir="${output_dir}"
  while [[ "${parent_dir}" != "/" ]]; do
    if [[ "${parent_dir##*/}" == "output" || "${parent_dir##*/}" == "figs" ]]; then
      output_dir="$(dirname -- "${parent_dir}")"
      break
    fi
    parent_dir="$(dirname -- "${parent_dir}")"
  done
  output_dir="${output_dir%/}/output"
  BRTP_OUT="${output_dir}/${base}_brtp.txt"
fi
mkdir -p -- "$(dirname -- "${BRTP_OUT}")"

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  # converter is absurdly tiny, but keep build behavior consistent with other wrappers
  "${ROOT_DIR}/build_gravmag_tools.sh"
fi

echo "+ step2 convert: ./gravmag_xyz_to_brtp ${XYZ_IN} ${BRTP_OUT}"
"${ROOT_DIR}/gravmag_xyz_to_brtp" "${XYZ_IN}" "${BRTP_OUT}"

echo "Done: ${BRTP_OUT}"
