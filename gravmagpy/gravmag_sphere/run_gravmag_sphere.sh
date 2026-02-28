#!/usr/bin/env bash
# ----------------------------------------------------------------------
# gravmag_sphere_brtp runner
#
# Usage:
#   ./run_gravmag_sphere.sh <R_sphere_km> <gauss_coeff_file> <input_in_file> [output_file]
# ----------------------------------------------------------------------

set -euo pipefail

# program names
SRC="gravmag_sphere_brtp.f"
EXE="gravmag_sphere_brtp"

# Args
R_SPHERE_KM="${1:-}"
GAUSS_FILE="${2:-}"
INPUT_FILE="${3:-}"
OUTFILE="${4:-}"

if [[ -z "${R_SPHERE_KM}" || -z "${GAUSS_FILE}" || -z "${INPUT_FILE}" ]]; then
  echo "ERROR: missing required arguments."
  echo "Usage: $0 <R_sphere_km> <gauss_coeff_file> <input_in_file> [output_file]"
  exit 2
fi

# sanity check for radius
if ! [[ "${R_SPHERE_KM}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: R_sphere_km must be a positive number (km). Got: ${R_SPHERE_KM}"
  exit 2
fi

# Files exist?
if [[ ! -f "${SRC}" ]]; then
  echo "ERROR: Fortran source not found: ${SRC}"
  echo "       (run from the directory containing ${SRC})"
  exit 2
fi

if [[ ! -f "${GAUSS_FILE}" ]]; then
  echo "ERROR: Gauss coefficient file not found: ${GAUSS_FILE}"
  exit 2
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "ERROR: Input .in file not found: ${INPUT_FILE}"
  exit 2
fi

# Clean staging
rm -f fort.1 fort.2 fort.4 sphout sphereiiscratch
rm -f "${EXE}" 2>/dev/null || true

# Compile
echo "+ gfortran -O2 -ffixed-form -std=legacy -ffixed-line-length-none -o ${EXE} ${SRC}"
gfortran -O2 -ffixed-form -std=legacy -ffixed-line-length-none \
  -o "${EXE}" "${SRC}"

# Stage unit files (Fortran opens unit 1 and 4 as fort.1 / fort.4)
echo "+ cp -f \"${GAUSS_FILE}\" fort.1"
cp -f "${GAUSS_FILE}" fort.1

echo "+ cp -f \"${INPUT_FILE}\" fort.4"
cp -f "${INPUT_FILE}" fort.4

# Run
echo "+ ./${EXE} ${R_SPHERE_KM}"
./"${EXE}" "${R_SPHERE_KM}"

# Collect outputs
echo "Done."
echo "Produced:"
ls -l fort.2 sphout sphereiiscratch 2>/dev/null || true

if [[ -n "${OUTFILE}" ]]; then
  if [[ -f fort.2 ]]; then
    echo "+ cp -f fort.2 \"${OUTFILE}\""
    cp -f fort.2 "${OUTFILE}"
    echo "Copied fort.2 -> ${OUTFILE}"
  else
    echo "WARNING: fort.2 not found (nfile may be 0 in input)."
  fi
fi

# clean up
rm fort.*
rm sphereiiscratch
rm sphout