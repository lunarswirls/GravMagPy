# SHTOOLS vs GravMagSphere (Spectral) Summary

## Scope

This compares:

- `shtools_lsq` (SHTOOLS least-squares spherical-harmonic fit)
- `fortran_spectral` (GravMagSphere Gauss/spectral solver)

using outputs already generated in this repo:

- `diagnostics/external_solver_comparison.csv`
- `diagnostics/gravmag_vs_shtools_residuals.csv`

All values below are from the current comparison setup (`lmax=24` for the external-comparison tables unless noted otherwise)

## Core method differences

### 1) How coefficients are estimated

- **SHTOOLS**:
  - Uses a spherical-harmonic least-squares expansion (`SHExpandLSQ`) from sampled data
  - Solves a linear inverse problem directly in coefficient space
- **GravMagSphere (spectral)**:
  - Builds/solves its own Gauss-coefficient system with explicit regularization controls (`reg_lambda`, `reg_power`) and geometry-specific preprocessing
  - Coefficients are then used to expand field components

### 2) Regularization and stability controls

- **SHTOOLS**:
  - Baseline LSQ behavior is relatively stable across test cases at `lmax=24`
- **GravMagSphere spectral**:
  - Accuracy is very sensitive to `lmax`, `reg_lambda`, and `reg_power`
  - Under-tuned settings can produce large residuals, especially for complex/sharp sources :(

### 3) Behavior near sharp boundaries

- Both approaches show Gibbs/ringing-type artifacts near discontinuities
- SHTOOLS edge artifacts appear as localized oscillatory bands around sharp source boundaries
- GravMagSphere shows stronger mismatch in difficult cases when coefficient fit/regularization is not tuned for that geometry

## Differences from current runs

From `diagnostics/external_solver_comparison.csv` (16 cases):

- Median RMSE(Btot):
  - `fortran_spectral`: **10.791**
  - `shtools_lsq`: **6.301**
- Case wins on RMSE(Btot):
  - `shtools_lsq` better in **12/16** cases
  - `fortran_spectral` better in **4/16** cases :(
- Median ratio `fortran_spectral / shtools_lsq` (RMSE(Btot)): **1.513**
  - GravMagSphere spectral error is ~1.5x SHTOOLS error at this setup :(

Heavy outliers in complex case (`gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge`):

- `fortran_spectral` RMSE(Btot): **10726.826**
- `shtools_lsq` RMSE(Btot): **278.470**

From direct GravMag-vs-SHTOOLS residuals (`diagnostics/gravmag_vs_shtools_residuals.csv`):

- Median RMSE(Btot) across cases: **7.807**
- Worst-case RMSE(Btot): **10713.753** (OOF)

Even with the same nominal `lmax`, not solving identical inverse problems...

- SHTOOLS LSQ minimizes SH fit residuals directly from sampled fields
- GravMagSphere spectral performs equivalent-source physics, joint-component fitting, and custom regularization before coefficient estimate

Case-dependent flip-flopping:

- Smooth/simple anomalies seem to favor direct LSQ efficiency and low residuals
- Sharp or highly structured boundaries can expose SH truncation limits in both models, but with different ringing signatures... Not sure why yet
- Multi-body results depend on how solver allocates shared low-degree coefficients across interacting sources
  - Stupid nonlinearity :(
