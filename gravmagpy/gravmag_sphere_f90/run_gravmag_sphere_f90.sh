#!/usr/bin/env bash
# ----------------------------------------------------------------------
# gravmag_sphere_brtp runner (modernized build, .mod in subdir)
#
# Usage:
#   ./run_gravmag_sphere.sh <R_sphere_km> <input_in_file> [output_file] [--debug]
#
# Debug options:
#   DEBUG=1 ./run_gravmag_sphere.sh ...
#   ./run_gravmag_sphere.sh ... --debug
# ----------------------------------------------------------------------

set -euo pipefail

# ----------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------
R_SPHERE_KM="${1:-}"
INPUT_FILE="${2:-}"
OUTFILE="${3:-}"
DEBUG_FLAG="${4:-}"

if [[ -z "${R_SPHERE_KM}" || -z "${INPUT_FILE}" ]]; then
  echo "ERROR: missing required arguments."
  echo "Usage: $0 <R_sphere_km> <input_in_file> [output_file] [--debug]"
  exit 2
fi

# Enable debug if env var or flag is present
DEBUG_MODE=0
if [[ "${DEBUG:-0}" == "1" || "${DEBUG_FLAG:-}" == "--debug" ]]; then
  DEBUG_MODE=1
fi

# ----------------------------------------------------------------------
# Sanity checks
# ----------------------------------------------------------------------
if ! [[ "${R_SPHERE_KM}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: R_sphere_km must be a positive number (km). Got: ${R_SPHERE_KM}"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: Input .in file not found: ${INPUT_FILE}"
  exit 2
fi

# Default output file if not provided
if [[ -z "${OUTFILE}" || "${OUTFILE}" == "--debug" ]]; then
  base="$(basename "${INPUT_FILE}")"
  base="${base%.*}"
  OUTFILE="${base}.txt"
fi

# Debug output name
BASE_INPUT="${INPUT_FILE%.in}"
DEBUG_OUT="${BASE_INPUT}.out"

# ----------------------------------------------------------------------
# Build directories
# ----------------------------------------------------------------------
MODDIR="mod"
mkdir -p "${MODDIR}"

# ----------------------------------------------------------------------
# Clean staging
# ----------------------------------------------------------------------
rm -f gauss.data input.in fort.2 sphout sphereiiscratch
rm -f gravmag_sphere_brtp *.o
rm -f "${MODDIR}"/*.mod 2>/dev/null || true

# ----------------------------------------------------------------------
# Compile
# ----------------------------------------------------------------------
echo "+ gfortran -std=f2008 -O2 -J${MODDIR} -I${MODDIR} \\"
echo "    gravmag_sphere_state.f90 \\"
echo "    gravmag_sphere_subs.f90 \\"
echo "    gravmag_sphere_physics.f90 \\"
echo "    gravmag_sphere_brtp.f90 \\"
echo "    -o gravmag_sphere_brtp"

gfortran -std=f2008 -O2 \
  -J"${MODDIR}" -I"${MODDIR}" \
  gravmag_sphere_state.f90 \
  gravmag_sphere_subs.f90 \
  gravmag_sphere_physics.f90 \
  gravmag_sphere_brtp.f90 \
  -o gravmag_sphere_brtp

# ----------------------------------------------------------------------
# Stage runtime files
# ----------------------------------------------------------------------
cp -f "${INPUT_FILE}" input.in

# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------
./gravmag_sphere_brtp "${R_SPHERE_KM}" input.in "${OUTFILE}"

# ----------------------------------------------------------------------
# Collect outputs
# ----------------------------------------------------------------------
echo "Done."
echo "Produced:"
ls -l "${OUTFILE}" fort.2 sphout sphereiiscratch 2>/dev/null || true

# Back-compat: fort.2 → OUTFILE if needed
if [[ ! -f "${OUTFILE}" && -f fort.2 ]]; then
  cp -f fort.2 "${OUTFILE}"
  echo "Copied fort.2 -> ${OUTFILE}"
fi

# ----------------------------------------------------------------------
# Debug handling
# ----------------------------------------------------------------------
if [[ "${DEBUG_MODE}" -eq 1 && -f sphout ]]; then
  cp -f sphout "${DEBUG_OUT}"
  echo "DEBUG file saved -> ${DEBUG_OUT}"
fi

# ----------------------------------------------------------------------
# Cleanup staging
# ----------------------------------------------------------------------
rm -f input.in fort.*
rm -f sphereiiscratch

# Only remove sphout if not debugging
if [[ "${DEBUG_MODE}" -eq 0 ]]; then
  rm -f sphout
fi