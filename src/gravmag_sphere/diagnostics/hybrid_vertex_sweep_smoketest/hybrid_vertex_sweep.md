# Hybrid Vertex Sweep

## Setup
- shapes: `regular,wavy`
- vertices: `4,8`
- spectral seeds: lmax `18`, reg_lambda `0.2`, reg_power `4.0`, auto_mode `0`
- hybrid params: band `1.5` deg, transition `0.75` deg, complex_vertex_threshold(test setting) `12`
- recommendation criteria: improvement >= `0.100`, runtime_ratio <= `2.500`, fraction_better >= `0.500`

## Recommendation
- recommended `complex_vertex_threshold`: **8**

## By Vertex (mean across shapes)

| n_vertices | mean RMSE improvement frac | mean runtime ratio (hyb/spec) | fraction hybrid better | fraction classified complex |
|---:|---:|---:|---:|---:|
| 4 | 0.0000 | 1.0100 | 0.0000 | 0.0000 |
| 8 | 0.0000 | 1.0011 | 0.0000 | 0.0000 |

## Per Case

| case | input_complex | hybrid_active | RMSEspec(Btot) | RMSEhyb(Btot) | improvement frac | runtime ratio |
|---|---:|---:|---:|---:|---:|---:|
| regular_v004 | 0 | 0 | 1.8607e+01 | 1.8607e+01 | 0.0000 | 1.0258 |
| regular_v008 | 0 | 0 | 2.0039e+01 | 2.0039e+01 | 0.0000 | 1.0004 |
| wavy_v004 | 0 | 0 | 8.2390e+01 | 8.2390e+01 | 0.0000 | 0.9941 |
| wavy_v008 | 0 | 0 | 7.9943e+01 | 7.9943e+01 | 0.0000 | 1.0019 |

## Files
- detailed csv: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep_smoketest/hybrid_vertex_sweep_detailed.csv`
- by-vertex csv: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep_smoketest/hybrid_vertex_sweep_by_vertex.csv`
- plot: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep_smoketest/hybrid_vertex_sweep.png`

## Methodology Note

This smoketest intentionally used a tiny vertex set (`4,8`) to validate pipeline mechanics, not to identify a stable production threshold.

- No hybrid gain here is expected because both shapes are classified simple under the configured complexity rules.
- Use the full sweep report for threshold selection and model-difference interpretation.
