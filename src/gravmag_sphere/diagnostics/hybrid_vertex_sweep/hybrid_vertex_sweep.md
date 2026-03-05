# Hybrid Vertex Sweep

## Setup
- shapes: `regular,wavy`
- vertices: `4,6,8,10,12,14,16,20,24,30,36,48`
- spectral seeds: lmax `18`, reg_lambda `0.2`, reg_power `4.0`, auto_mode `0`
- hybrid params: band `1.5` deg, transition `0.75` deg, complex_vertex_threshold(test setting) `12`
- recommendation criteria: improvement >= `0.100`, runtime_ratio <= `2.500`, fraction_better >= `0.500`

## Recommendation
- recommended `complex_vertex_threshold`: **12**

## By Vertex (mean across shapes)

| n_vertices | mean RMSE improvement frac | mean runtime ratio (hyb/spec) | fraction hybrid better | fraction classified complex |
|---:|---:|---:|---:|---:|
| 4 | 0.0000 | 0.9833 | 0.0000 | 0.0000 |
| 6 | 0.0000 | 1.0019 | 0.0000 | 0.0000 |
| 8 | 0.0000 | 0.9924 | 0.0000 | 0.0000 |
| 10 | 0.0278 | 1.0277 | 0.0000 | 0.5000 |
| 12 | 0.1570 | 0.9996 | 0.5000 | 0.5000 |
| 14 | 0.1323 | 1.0161 | 0.5000 | 1.0000 |
| 16 | 0.2980 | 1.0029 | 1.0000 | 1.0000 |
| 20 | 0.1089 | 1.0138 | 0.5000 | 1.0000 |
| 24 | 0.2463 | 1.0006 | 1.0000 | 1.0000 |
| 30 | 0.1260 | 1.0222 | 0.5000 | 1.0000 |
| 36 | 0.2336 | 1.0086 | 1.0000 | 1.0000 |
| 48 | 0.1680 | 1.0146 | 0.5000 | 1.0000 |

## Per Case

| case | input_complex | hybrid_active | RMSEspec(Btot) | RMSEhyb(Btot) | improvement frac | runtime ratio |
|---|---:|---:|---:|---:|---:|---:|
| regular_v004 | 0 | 0 | 1.8607e+01 | 1.8607e+01 | 0.0000 | 0.9818 |
| regular_v006 | 0 | 0 | 1.6103e+01 | 1.6103e+01 | 0.0000 | 1.0036 |
| regular_v008 | 0 | 0 | 2.0039e+01 | 2.0039e+01 | 0.0000 | 1.0020 |
| regular_v010 | 0 | 0 | 1.9013e+01 | 1.9013e+01 | 0.0000 | 1.0350 |
| regular_v012 | 0 | 0 | 1.9172e+01 | 1.9172e+01 | 0.0000 | 0.9812 |
| regular_v014 | 1 | 1 | 1.8382e+01 | 1.3358e+01 | 0.2733 | 1.0260 |
| regular_v016 | 1 | 1 | 1.8910e+01 | 1.4106e+01 | 0.2541 | 0.9909 |
| regular_v020 | 1 | 1 | 1.8951e+01 | 1.4074e+01 | 0.2573 | 1.0328 |
| regular_v024 | 1 | 1 | 1.9428e+01 | 1.3833e+01 | 0.2880 | 0.9938 |
| regular_v030 | 1 | 1 | 1.9082e+01 | 1.3633e+01 | 0.2856 | 1.0466 |
| regular_v036 | 1 | 1 | 1.9474e+01 | 1.3828e+01 | 0.2899 | 1.0165 |
| regular_v048 | 1 | 1 | 1.9474e+01 | 1.3845e+01 | 0.2891 | 1.0092 |
| wavy_v004 | 0 | 0 | 8.2390e+01 | 8.2390e+01 | 0.0000 | 0.9848 |
| wavy_v006 | 0 | 0 | 8.6739e+01 | 8.6739e+01 | 0.0000 | 1.0003 |
| wavy_v008 | 0 | 0 | 7.9943e+01 | 7.9943e+01 | 0.0000 | 0.9828 |
| wavy_v010 | 1 | 1 | 4.3200e+01 | 4.0798e+01 | 0.0556 | 1.0203 |
| wavy_v012 | 1 | 1 | 1.6978e+01 | 1.1646e+01 | 0.3141 | 1.0179 |
| wavy_v014 | 1 | 1 | 2.5612e+01 | 2.5836e+01 | -0.0088 | 1.0061 |
| wavy_v016 | 1 | 1 | 5.2544e+01 | 3.4575e+01 | 0.3420 | 1.0149 |
| wavy_v020 | 1 | 1 | 4.0430e+01 | 4.2030e+01 | -0.0396 | 0.9947 |
| wavy_v024 | 1 | 1 | 1.7411e+01 | 1.3848e+01 | 0.2046 | 1.0075 |
| wavy_v030 | 1 | 1 | 4.1181e+01 | 4.2561e+01 | -0.0335 | 0.9978 |
| wavy_v036 | 1 | 1 | 1.8273e+01 | 1.5035e+01 | 0.1772 | 1.0008 |
| wavy_v048 | 1 | 1 | 1.8169e+01 | 1.7316e+01 | 0.0469 | 1.0201 |

## Files
- detailed csv: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep_detailed.csv`
- by-vertex csv: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep_by_vertex.csv`
- plot: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/hybrid_vertex_sweep/hybrid_vertex_sweep.png`

## Methodology Note: Why Hybrid Helps Above the Threshold

- Spectral-only mode represents fields with finite SH bandwidth, which smooths discontinuities and can ring near polygon boundaries.
- Hybrid mode replaces spectral prediction with direct kernel evaluation inside a boundary band and blends back outside that band.

Why the recommendation lands near 12 vertices in this sweep:

- Below the threshold, shapes are simple enough that spectral smoothing cost is low and hybrid offers little gain.
- Above the threshold, geometric complexity and boundary curvature increase high-degree demand; local direct correction recovers edge behavior with minor runtime cost.
- Improvement is not strictly monotonic because waviness and boundary orientation change how much unresolved energy sits in high degrees.
