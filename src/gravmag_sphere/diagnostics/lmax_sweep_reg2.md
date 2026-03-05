# Spherical Harmonic lmax Sweep

- Compiler: `gfortran`
- Example: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in`
- Edge band: `0.150` deg
- Far band: `0.400` deg
- Ringing thresholds: edge_amp <= `2.0`, edge_overshoot <= `0.15`

| lmax | Wall (s) | RMSE | Edge amp | Edge overshoot | Ringing OK |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.337 | 4.116e+00 | 92.681 | 48.573 | no |
| 16 | 0.705 | 4.115e+00 | 93.423 | 48.567 | no |
| 24 | 3.330 | 4.780e+00 | 9.873 | 50.842 | no |

Recommended max lmax: none (thresholds too strict for this setup).
Caution: this is empirical for the selected case, solver regularization, mesh, and thresholds.

## Why Solutions Change with `lmax` and Alternate Regularization

- This run uses a different regularization setting, so high-degree damping differs from the default-fit sweep.
- As `lmax` increases, added degrees can improve edge localization but can also reintroduce oscillatory boundary artifacts.
- Differences relative to the default sweep come from how the changed regularization redistributes energy across degree.
- The key takeaway is still joint tuning: `lmax` and damping must be selected together for stable edge behavior.
