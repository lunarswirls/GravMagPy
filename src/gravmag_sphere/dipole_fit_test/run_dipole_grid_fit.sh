#!/usr/bin/env bash
set -euo pipefail

# Dipole-grid inversion runner
# Usage:
#   ./dipole_fit_test/run_dipole_grid_fit.sh <R_sphere_km> <obs_csvs> [predictions_csv] [dipoles_csv] \
#     [source_depth_km] [source_dlat_deg] [source_dlon_deg] [reg_lambda] [source_nr] \
#     [source_dr_km] [lat_pad_deg] [lon_pad_deg] [max_memory_mib]
#
# <obs_csvs> is a comma-separated list of CSV files. Each CSV must contain:
#   - lon/lat columns
#   - altitude_m|altitude_km|radius_m|radius_km
#   - Bx/By/Bz in either nT or Tesla

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
R_SPHERE_KM="${1:-}"
OBS_CSVS="${2:-}"
PRED_CSV="${3:-}"
DIPOLE_CSV="${4:-}"
SOURCE_DEPTH_KM="${5:-0}"
SOURCE_DLAT_DEG="${6:-0}"
SOURCE_DLON_DEG="${7:-0}"
REG_LAMBDA="${8:-1.0e-4}"
SOURCE_NR="${9:-1}"
SOURCE_DR_KM="${10:-0}"
LAT_PAD_DEG="${11:--1}"
LON_PAD_DEG="${12:--1}"
MAX_MEMORY_MIB="${13:-768}"

if [[ -z "${R_SPHERE_KM}" || -z "${OBS_CSVS}" ]]; then
  echo "Usage: $0 <R_sphere_km> <obs_csvs> [predictions_csv] [dipoles_csv] [source_depth_km] [source_dlat_deg] [source_dlon_deg] [reg_lambda] [source_nr] [source_dr_km] [lat_pad_deg] [lon_pad_deg] [max_memory_mib]"
  exit 2
fi

IFS=',' read -r FIRST_CSV _ <<< "${OBS_CSVS}"
if [[ -z "${FIRST_CSV}" || ! -f "${FIRST_CSV}" ]]; then
  echo "ERROR: first observation CSV not found: ${FIRST_CSV}"
  exit 2
fi

if [[ -z "${PRED_CSV}" || -z "${DIPOLE_CSV}" ]]; then
  base="$(basename "${FIRST_CSV}")"
  base="${base%.csv}"
  mkdir -p "${ROOT_DIR}/output"
  [[ -z "${PRED_CSV}" ]] && PRED_CSV="${ROOT_DIR}/output/${base}_dipole_fit_predictions.csv"
  [[ -z "${DIPOLE_CSV}" ]] && DIPOLE_CSV="${ROOT_DIR}/output/${base}_dipole_fit_dipoles.csv"
fi

"${ROOT_DIR}/build_gravmag_tools.sh"

echo "+ running dipole-grid inversion"
echo "+ ./dipole_fit_test/gravmag_sphere_dipole_grid_fit ${R_SPHERE_KM} ${OBS_CSVS} ${PRED_CSV} ${DIPOLE_CSV} ${SOURCE_DEPTH_KM} ${SOURCE_DLAT_DEG} ${SOURCE_DLON_DEG} ${REG_LAMBDA} ${SOURCE_NR} ${SOURCE_DR_KM} ${LAT_PAD_DEG} ${LON_PAD_DEG} ${MAX_MEMORY_MIB}"
"${SCRIPT_DIR}/gravmag_sphere_dipole_grid_fit" \
  "${R_SPHERE_KM}" "${OBS_CSVS}" "${PRED_CSV}" "${DIPOLE_CSV}" \
  "${SOURCE_DEPTH_KM}" "${SOURCE_DLAT_DEG}" "${SOURCE_DLON_DEG}" "${REG_LAMBDA}" \
  "${SOURCE_NR}" "${SOURCE_DR_KM}" "${LAT_PAD_DEG}" "${LON_PAD_DEG}" "${MAX_MEMORY_MIB}"

echo "Prediction CSV: ${PRED_CSV}"
echo "Dipole CSV: ${DIPOLE_CSV}"
