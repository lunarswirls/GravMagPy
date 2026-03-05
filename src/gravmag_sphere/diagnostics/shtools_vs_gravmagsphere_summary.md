# SHTOOLS vs GravMagSphere (Fortran Spectral) Summary

## Scope of this summary

This compares:

- `shtools_lsq` (SHTOOLS least-squares spherical-harmonic fit)
- `fortran_spectral` (GravMagSphere Gauss/spectral solver)

using outputs already generated in this repo:

- `diagnostics/external_solver_comparison.csv`
- `diagnostics/gravmag_vs_shtools_residuals.csv`

All values below are from the current comparison setup (`lmax=24` for the external-comparison tables unless noted otherwise).

## Core method differences

### 1) How coefficients are estimated

- **SHTOOLS**:
  - Uses a spherical-harmonic least-squares expansion (`SHExpandLSQ`) from sampled data.
  - Solves a linear inverse problem directly in coefficient space.
- **GravMagSphere spectral (Fortran)**:
  - Builds/solves its own Gauss-coefficient system with explicit regularization controls (`reg_lambda`, `reg_power`) and geometry-specific preprocessing in the GravMag pipeline.
  - Coefficients are then used to synthesize field components.

### 2) Regularization and stability controls

- **SHTOOLS**:
  - Baseline LSQ behavior in this workflow is relatively stable across many test cases at `lmax=24`.
  - Does not use the GravMagSphere-specific `reg_lambda/reg_power` controls.
- **GravMagSphere spectral**:
  - Accuracy is very sensitive to `lmax`, `reg_lambda`, and `reg_power`.
  - Under-tuned settings can produce large residuals, especially for complex/sharp sources.

### 3) Behavior near sharp boundaries

- Both approaches can show Gibbs/ringing-type artifacts near discontinuities.
- In your residual plots, SHTOOLS-related edge artifacts appear as localized oscillatory bands around sharp source limits.
- GravMagSphere can show stronger mismatch in difficult cases when coefficient fit/regularization is not tuned for that geometry.

## Quantitative differences from current runs

From `diagnostics/external_solver_comparison.csv` (16 cases):

- Median RMSE(Btot):
  - `fortran_spectral`: **10.791**
  - `shtools_lsq`: **6.301**
- Case wins on RMSE(Btot):
  - `shtools_lsq` better in **12/16** cases
  - `fortran_spectral` better in **4/16** cases
- Median ratio `fortran_spectral / shtools_lsq` (RMSE(Btot)): **1.513**
  - Typical case: Fortran spectral error is ~1.5x SHTOOLS error at this setup.

Heavy outliers in complex case (`gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge`):

- `fortran_spectral` RMSE(Btot): **10726.826**
- `shtools_lsq` RMSE(Btot): **278.470**

From direct GravMag-vs-SHTOOLS residuals (`diagnostics/gravmag_vs_shtools_residuals.csv`):

- Median RMSE(Btot) across cases: **7.807**
- Worst-case RMSE(Btot): **10713.753** (same complexlarge family)

## Pros and cons

### SHTOOLS (`shtools_lsq`)

#### Pros

- Strong baseline robustness in your current test set.
- Better RMSE(Btot) in most cases at `lmax=24`.
- Mature spherical-harmonic tooling and predictable LSQ behavior.

#### Cons

- Still exhibits edge/ringing artifacts around sharp discontinuities.
- As a pure SH model, discontinuous source geometries need higher `lmax` or filtering/tapering to suppress ringing.
- May require extra post-processing choices to control edge behavior for production maps.

### GravMagSphere spectral (`fortran_spectral`)

#### Pros

- Native integration with your GravMagSphere workflow and formats.
- Exposes explicit regularization controls (`reg_lambda`, `reg_power`) for targeted tuning.
- Can be competitive or better in a subset of cases when tuned appropriately.

#### Cons

- More sensitive to parameter choice; can fail badly on complex/sharp geometries if under-tuned.
- Current default-like settings used in comparison show larger residuals in many cases.
- Requires wider parameter sweeps (as you are doing) to identify stable operating regions.

## Practical takeaway

- If you need a stable baseline now, **SHTOOLS LSQ is currently the safer default** in this dataset.
- If you need GravMagSphere-native production runs, **use the operating-map sweep results** to choose `(lmax, reg_lambda, reg_power)` by runtime/accuracy target rather than fixed defaults.
- For sharp-edged sources, treat both solvers as spectral approximations and plan for anti-ringing strategy (higher `lmax`, localized smoothing/taper, or hybrid validation against direct solver).

## Source files

- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/external_solver_comparison.md`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/external_solver_comparison.csv`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/gravmag_vs_shtools_residuals.md`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/gravmag_vs_shtools_residuals.csv`

## Additional Explanation: Why Case-by-Case Rankings Flip

Even with the same nominal `lmax`, the models are not solving identical inverse problems.

- SHTOOLS LSQ minimizes SH fit residuals directly from sampled fields.
- GravMagSphere spectral carries equivalent-source physics, joint-component fitting, and custom regularization into the coefficient estimate.

This produces case-dependent tradeoffs:

- Smooth/simple anomalies often favor direct SH LSQ efficiency and low residuals.
- Sharp or highly structured boundaries can expose SH truncation limits in both models, but with different ringing signatures.
- Multi-body scenes depend on how each solver allocates shared low-degree coefficients across interacting sources.
- Local edge corrections/hybrid blending in GravMagSphere can reduce boundary misfit without making the global SH model identical to SHTOOLS.
