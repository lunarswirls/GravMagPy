# How `external_solver_residuals` Are Computed

This note documents exactly how residuals in
`src/gravmag_sphere/diagnostics/external_solver_residuals` are generated.

Implementation source: `src/gravmag_sphere/external_solver_compare.py`.

## 1) Reference Baseline

For each example case, the reference field is the **direct solver** output:

1. Run:
   - `./run_input_to_xyz.sh direct ...`
2. Convert:
   - `./run_xyz_to_brtp.sh ...`
3. Load and aggregate by `(lon, lat)` using `load_brtp_aggregate(...)`.

Residuals are always computed against this direct baseline.

## 2) Predicted Solver Field

For each case, three predictors are compared to the direct baseline:

1. `fortran_spectral`
   - Runs `./run_input_to_xyz.sh spectral ...` with configured `lmax`, `reg_lambda`, `reg_power`.
   - Converted to BRTP and aggregated by `(lon, lat)`.
2. `scipy_sh_lsq`
   - Builds a real SH design matrix and solves least squares (optionally ridge-regularized).
3. `shtools_lsq`
   - Uses `pyshtools.expand.SHExpandLSQ` + `pyshtools.expand.LSQ_G`.

For all predictors, `Btot` is recomputed as:

`Btot = sqrt(Br^2 + Btheta^2 + Bphi^2)`

## 3) Grid Alignment Before Residuals

Baseline and prediction are aligned by the **intersection** of rounded coordinates:

- key: `(round(lon, 6), round(lat, 6))`
- overlap only: samples present in both baseline and prediction are used.

This is handled by `align_component_maps(...)`.

## 4) Residual Definition

Per component (`Br`, `Btheta`, `Bphi`, `Btot`), residual is:

`residual = prediction - direct_baseline`

So positive values mean the tested solver predicts larger values than direct baseline.

## 5) Reported Metrics

Computed per component from residual vector `r = y_pred - y_true`:

1. `RMSE = sqrt(mean(r^2))`
2. `MAE = mean(abs(r))`
3. `max_abs = max(abs(r))`
4. `R2 = 1 - sum(r^2) / sum((y_true - mean(y_true))^2)`

If `sum((y_true - mean(y_true))^2) <= 0`, `R2` is reported as `NaN`.

## 6) Residual Plot Construction

Each residual PNG in `external_solver_residuals` contains four panels:

1. `Br residual`
2. `Btheta residual`
3. `Bphi residual`
4. `Btot residual`

Plot details:

- residual field is gridded with `to_grid(...)`
- colormap: `RdBu_r`
- color limits are symmetric per panel:
  - `vmin = -max(abs(residual))`
  - `vmax = +max(abs(residual))`
- colorbar label: `Pred - baseline`

## 7) Regeneration Command

From `src/gravmag_sphere`:

```bash
.venv_compare/bin/python external_solver_compare.py compare \
  --formatted-root external_solver_inputs \
  --report-md diagnostics/external_solver_comparison.md \
  --report-csv diagnostics/external_solver_comparison.csv \
  --residual-dir diagnostics/external_solver_residuals \
  --lmax 24 \
  --ridge-lambda 1e-8 \
  --rsphere-km 1737.4 \
  --refine-factor 2 \
  --ntheta-fit 72 \
  --nphi-fit 144 \
  --reg-lambda 0.2 \
  --reg-power 4.0
```


## 7) Methodological Interpretation of Model-to-Model Differences

These residuals are not only numerical error; they also encode model-assumption mismatch.

- GravMagSphere direct uses physically discretized bodies and kernel summation.
- GravMagSphere spectral projects that physics into a truncated SH basis with regularization.
- SciPy/SHTOOLS LSQ rows in this workflow solve a data-fit SH inverse model directly.

Typical pattern interpretation:

- Narrow red/blue edge bands: Gibbs-type behavior from finite SH bandwidth at discontinuous boundaries.
- Broad smooth bias over an anomaly: low-degree dominance from regularization or underfit `lmax`.
- Strong `Bphi` mismatch near source corners: azimuthal component sensitivity to basis phase and component-wise inversion.
- Multi-body misfit halos: coefficient competition between nearby sources and imperfect separation in low-order harmonics.

Because these methods do not share identical priors, equal `lmax` does not imply equivalent solutions.
