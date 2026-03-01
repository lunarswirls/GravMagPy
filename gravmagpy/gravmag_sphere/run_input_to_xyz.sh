#!/usr/bin/env bash
set -euo pipefail

# Step 1: input .in -> solver XYZ output.
#
# Usage:
#   ./run_input_to_xyz.sh <solver> <R_sphere_km> <input_in_file> \
#     [xyz_output_file] [refine_factor] [lmax] [ntheta_fit] [nphi_fit] \
#     [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr]
#
# solver:
#   - spectral : gravmag_sphere_gauss
#   - direct   : gravmag_sphere_bxyz
#
# output:
#   - writes one XYZ table (body_id lon lat Fx Fy Fz Ftot)
#   - does not perform spherical component conversion.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SOLVER="${1:-spectral}"
R_SPHERE_KM="${2:-}"
INPUT_FILE="${3:-}"
XYZ_OUT="${4:-}"
REFINE_FACTOR="${5:-2}"
LMAX="${6:-24}"
NTHETA_FIT="${7:-72}"
NPHI_FIT="${8:-144}"
REG_LAMBDA="${9:-0.2}"
REG_POWER="${10:-4.0}"
SOURCE_NLAT="${11:-0}"
SOURCE_NLON="${12:-0}"
SOURCE_NR="${13:-0}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "Usage: $0 <solver> <R_sphere_km> <input_in_file> [xyz_output_file] [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr]"
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

if [[ "${SKIP_BUILD:-0}" != "1" ]]; then
  # build once unless caller explicitly skipped it.
  "${ROOT_DIR}/build_gravmag_tools.sh"
fi

case "${SOLVER}" in
  spectral)
    # spectral pipeline: source mesh -> harmonic fit -> XYZ field table.
    echo "+ step1 spectral: ./gravmag_sphere_gauss ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${LMAX} ${REFINE_FACTOR} ${NTHETA_FIT} ${NPHI_FIT} ${REG_LAMBDA} ${REG_POWER} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
    "${ROOT_DIR}/gravmag_sphere_gauss" \
      "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" \
      "${LMAX}" "${REFINE_FACTOR}" "${NTHETA_FIT}" "${NPHI_FIT}" \
      "${REG_LAMBDA}" "${REG_POWER}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"
    ;;
  direct)
    # direct pipeline: source mesh -> direct Green's function summation.
    echo "+ step1 direct: ./gravmag_sphere_bxyz ${R_SPHERE_KM} ${INPUT_FILE} ${XYZ_OUT} ${REFINE_FACTOR} ${SOURCE_NLAT} ${SOURCE_NLON} ${SOURCE_NR}"
    "${ROOT_DIR}/gravmag_sphere_brtp" \
      "${R_SPHERE_KM}" "${INPUT_FILE}" "${XYZ_OUT}" \
      "${REFINE_FACTOR}" "${SOURCE_NLAT}" "${SOURCE_NLON}" "${SOURCE_NR}"
    ;;
  *)
    echo "ERROR: solver must be 'spectral' or 'direct' (got: ${SOLVER})"
    exit 2
    ;;
esac

echo "Done: ${XYZ_OUT}"
