#!/usr/bin/env bash
set -euo pipefail

# Direct solver runner (XYZ output).
# Usage:
#   ./run_gravmag_sphere_f90.sh <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]
#
# Equivalent executable call:
#   ./gravmag_sphere_brtp R input.in output_xyz.txt refine source_nlat source_nlon source_nr

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SPHERE_KM="${1:-}"
INPUT_FILE="${2:-}"
XYZ_OUT="${3:-}"
REFINE_FACTOR="${4:-2}"
SOURCE_NLAT="${5:-0}"
SOURCE_NLON="${6:-0}"
SOURCE_NR="${7:-0}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [source_nlat] [source_nlon] [source_nr]"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: input file not found: ${INPUT_FILE}"
  exit 2
fi

if [[ -z "${XYZ_OUT}" ]]; then
  base="$(basename "${INPUT_FILE}")"
  base="${base%.in}"
  mkdir -p "${ROOT_DIR}/output"
  XYZ_OUT="${ROOT_DIR}/output/${base}_xyz.txt"
fi

"${ROOT_DIR}/build_gravmag_tools.sh"

echo "+ running direct solver"
echo "+ ./gravmag_sphere_brtp ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${REFINE_FACTOR} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
"${ROOT_DIR}/gravmag_sphere_brtp" "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" "${REFINE_FACTOR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"

echo "Done: ${XYZ_OUT}"
