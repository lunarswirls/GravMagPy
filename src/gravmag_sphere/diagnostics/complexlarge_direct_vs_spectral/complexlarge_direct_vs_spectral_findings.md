# Complexlarge Direct vs Spectral Boundary Analysis

- run directory: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/complexlarge_direct_vs_spectral`
- side-by-side fields: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/complexlarge_direct_vs_spectral/complexlarge_direct_vs_spectral_side_by_side.png`
- residuals (spectral - direct): `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/complexlarge_direct_vs_spectral/complexlarge_spectral_minus_direct_residuals.png`
- residuals (spectral no-local-correction - direct): `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/complexlarge_direct_vs_spectral/complexlarge_spectral_nolocal_minus_direct_residuals.png`

## Metrics: Spectral - Direct

| Component | RMSE | MAE | MaxAbs |
|---|---:|---:|---:|
| bx | 1.661097e+02 | 5.165958e+01 | 1.461458e+04 |
| by | 1.545729e+02 | 4.387235e+01 | 9.867355e+03 |
| bz | 2.335184e+02 | 4.557321e+01 | 2.530574e+04 |
| btot | 3.189879e+02 | 8.936603e+01 | 2.974910e+04 |

## Metrics: Spectral(no local correction) - Direct

| Component | RMSE | MAE | MaxAbs |
|---|---:|---:|---:|
| bx | 1.376958e+02 | 4.966661e+01 | 1.459403e+04 |
| by | 1.334884e+02 | 4.179716e+01 | 9.867307e+03 |
| bz | 2.261194e+02 | 4.316234e+01 | 2.529955e+04 |
| btot | 2.919912e+02 | 8.446118e+01 | 2.972394e+04 |

## Metrics: Spectral(local) - Spectral(no local correction)

| Component | RMSE | MAE | MaxAbs |
|---|---:|---:|---:|
| bx | 9.244795e+01 | 6.165520e+00 | 5.874053e+03 |
| by | 8.513715e+01 | 5.130115e+00 | 5.897294e+03 |
| bz | 7.984852e+01 | 3.892112e+00 | 5.522899e+03 |
| btot | 1.458752e+02 | 8.023746e+00 | 8.780817e+03 |

## Edge-Region Contrast (using top-10% direct Btot gradient as boundary proxy)

| Model | MeanAbsResidual Near Edge | MeanAbsResidual Far From Edge | Near/Far Ratio |
|---|---:|---:|---:|
| spectral | 2.218059e+02 | 7.464443e+01 | 2.972 |
| spectral_nolocal | 1.727446e+02 | 7.464788e+01 | 2.314 |

## Interpretation

The visible boundary is primarily a high-gradient edge effect from approximating a sharp polygonal source with a finite-order spectral model. Residuals concentrate near edge-like zones much more strongly than away from edges.
Local edge correction changes the solution significantly near boundaries (see local-vs-nolocal metrics), which can make boundary contrast more visible even when broad ringing is reduced.

## Methodology Difference Behind the Visible Boundary

- Direct model uses physical source-kernel summation and preserves discontinuous boundary behavior.
- Spectral model uses finite SH representation, which cannot reproduce abrupt polygon edges without oscillatory spillover.
- Local edge correction modifies near-boundary amplitudes but cannot fully remove global truncation effects from the SH backbone.

Why boundaries remain prominent:

- High gradient at polygon limits injects high-degree spectral content beyond chosen bandwidth.
- Residual energy concentrates where boundary curvature changes rapidly.
- In this case, edge-zone mismatch dominates far-field mismatch, so the boundary appears visually emphasized in residual maps.
