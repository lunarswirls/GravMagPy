#!/usr/bin/env bash
set -euo pipefail

# Build all Fortran executables used in this workspace:
#   - gravmag_sphere_bxyz      (direct volume solver, XYZ output)
#   - gravmag_sphere_gauss     (spectral/gauss alternative)
#   - gravmag_xyz_to_brtp      (XYZ -> Br/Btheta/Bphi converter)
#
# Notes:
#   - this script is intentionally explicit (no makefile dependency),
#       so every compile command is visible and easy to edit :)
#   - module files are written to ./mod

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODDIR="${ROOT_DIR}/mod"
mkdir -p "${MODDIR}"

echo "+ building gravmag_sphere_bxyz"
# direct solver depends on three helper modules
gfortran -std=f2008 -O2 -fopenmp -J"${MODDIR}" -I"${MODDIR}" \
  "${ROOT_DIR}/gravmag_sphere_state.f90" \
  "${ROOT_DIR}/gravmag_sphere_subs.f90" \
  "${ROOT_DIR}/gravmag_sphere_physics.f90" \
  "${ROOT_DIR}/gravmag_sphere_bxyz.f90" \
  -o "${ROOT_DIR}/gravmag_sphere_bxyz"

echo "+ building gravmag_sphere_gauss"
# spectral solver is currently self-contained in one source file
gfortran -std=f2008 -O2 \
  "${ROOT_DIR}/gravmag_sphere_gauss.f90" \
  -o "${ROOT_DIR}/gravmag_sphere_gauss"

echo "+ building gravmag_xyz_to_brtp"
# xyz->brtp converter is also self-contained
gfortran -std=f2008 -O2 \
  "${ROOT_DIR}/gravmag_xyz_to_brtp.f90" \
  -o "${ROOT_DIR}/gravmag_xyz_to_brtp"

echo "Build complete."
