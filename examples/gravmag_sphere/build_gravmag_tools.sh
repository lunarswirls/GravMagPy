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
source_dir="${ROOT_DIR}/../../fortran"
mkdir -p "${MODDIR}"

echo "+ building gravmag_sphere_bxyz"
# direct solver depends on three helper modules
gfortran -std=f2008 -O2 -fopenmp -J"${MODDIR}" -I"${MODDIR}" \
  "${source_dir}/gravmag_paths.f90" \
  "${source_dir}/gravmag_sphere_state.f90" \
  "${source_dir}/gravmag_sphere_subs.f90" \
  "${source_dir}/gravmag_sphere_physics.f90" \
  "${source_dir}/gravmag_sphere_bxyz.f90" \
  -o "${ROOT_DIR}/gravmag_sphere_bxyz"

echo "+ building gravmag_sphere_gauss"
# share only path handling; numerical solvers remain independent
gfortran -std=f2008 -O2 -J"${MODDIR}" -I"${MODDIR}" \
  "${source_dir}/gravmag_paths.f90" \
  "${source_dir}/gravmag_sphere_gauss.f90" \
  -o "${ROOT_DIR}/gravmag_sphere_gauss"

echo "+ building gravmag_xyz_to_brtp"
# build the converter with the shared output path utility
gfortran -std=f2008 -O2 -J"${MODDIR}" -I"${MODDIR}" \
  "${source_dir}/gravmag_paths.f90" \
  "${source_dir}/gravmag_xyz_to_brtp.f90" \
  -o "${ROOT_DIR}/gravmag_xyz_to_brtp"

# compile the existing optional dipole-grid workflow too
gfortran -std=f2008 -O2 -J"${MODDIR}" -I"${MODDIR}" \
  "${source_dir}/gravmag_paths.f90" "${source_dir}/gravmag_sphere_dipole_grid_fit.f90" \
  -o "${ROOT_DIR}/dipole_fit_test/gravmag_sphere_dipole_grid_fit"

echo "Build complete."
