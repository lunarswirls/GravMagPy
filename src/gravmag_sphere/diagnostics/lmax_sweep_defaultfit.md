# Spherical Harmonic lmax Sweep

- Compiler: `gfortran`
- Example: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in`
- Edge band: `0.150` deg
- Far band: `0.400` deg
- Ringing thresholds: edge_amp <= `2.0`, edge_overshoot <= `0.15`

| lmax | Wall (s) | RMSE | Edge amp | Edge overshoot | Ringing OK |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.361 | 4.116e+00 | 92.747 | 48.573 | no |
| 16 | 1.127 | 4.114e+00 | 93.615 | 48.563 | no |
| 24 | 3.435 | 4.710e+00 | 11.087 | 49.792 | no |

Recommended max lmax: none (thresholds too strict for this setup).
Caution: this is empirical for the selected case, solver regularization, mesh, and thresholds.

## Why Solutions Change with `lmax`

- Increasing `lmax` expands the harmonic bandwidth and can represent sharper features.
- The same change also amplifies sensitivity to discontinuities, which can increase edge ringing when regularization is insufficient.
- Therefore error behavior is non-monotonic: some boundary metrics improve while overshoot can worsen.
- This sweep shows that `lmax` tuning must be paired with regularization/hybrid-edge strategy, not increased in isolation.
