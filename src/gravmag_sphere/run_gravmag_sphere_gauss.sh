#!/usr/bin/env bash
set -euo pipefail

# Gauss/spectral solver runner.
# Usage:
#   ./run_gravmag_sphere_gauss.sh <R_sphere_km> <input_in_file> [output_file] [lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr]
#
# This wrapper exposes all main spectral controls:
#   - lmax controls harmonic resolution,
#   - reg_lambda/reg_power control high-degree damping.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SPHERE_KM="${1:-}"
INPUT_FILE="${2:-}"
OUTFILE="${3:-}"
LMAX="${4:-24}"
REFINE="${5:-2}"
NTHETA_FIT="${6:-72}"
NPHI_FIT="${7:-144}"
REG_LAMBDA="${8:-0.2}"
REG_POWER="${9:-4.0}"
SOURCE_NLAT="${10:-0}"
SOURCE_NLON="${11:-0}"
SOURCE_NR="${12:-0}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <R_sphere_km> <input_in_file> [output_file] [lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr]"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: input file not found: ${INPUT_FILE}"
  exit 2
fi

if [[ -z "${OUTFILE}" ]]; then
  base="$(basename "${INPUT_FILE}")"
  base="${base%.in}"
  mkdir -p "${ROOT_DIR}/output"
  OUTFILE="${ROOT_DIR}/output/${base}_gauss.txt"
fi

"${ROOT_DIR}/build_gravmag_tools.sh"

echo "+ running gauss solver"
echo "+ ./gravmag_sphere_gauss ${R_SPHERE_KM} ${INPUT_FILE} ${OUTFILE} ${LMAX} ${REFINE} ${NTHETA_FIT} ${NPHI_FIT} ${REG_LAMBDA} ${REG_POWER} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
"${ROOT_DIR}/gravmag_sphere_gauss" "${R_SPHERE_KM}" "${INPUT_FILE}" "${OUTFILE}" "${LMAX}" "${REFINE}" "${NTHETA_FIT}" "${NPHI_FIT}" "${REG_LAMBDA}" "${REG_POWER}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"

echo "Done: ${OUTFILE}"
