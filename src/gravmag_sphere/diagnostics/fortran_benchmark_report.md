# Fortran Compiler Benchmark

- Example: `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in`

| Compiler | Case | Median wall (s) | Slowest stage | Largest memory bucket |
|---|---:|---:|---|---|
| gfortran | direct_polygon | 0.105 | field_accum_s (0.090s) | total_est_mib (0.31 MiB) |
| gfortran | spectral_polygon | 3.430 | fit_assembly_s (2.790s) | total_est_mib (3.58 MiB) |

Notes: stage and memory values come from DIAG instrumentation (`GRAVMAG_DIAGNOSTICS=1`).

## Performance

- `direct_polygon` scales mainly with number of source elements times observation cells (kernel summation dominated)
- `spectral_polygon` adds fit-grid assembly and linear solve costs that scale strongly with harmonic degree and sample count
- This is why spectral can use more memory and spend most time in fit assembly/solve even when field evaluation itself is fast
- Solution differences between direct and spectral are driven by basis truncation and regularization, while runtime differences are driven by inversion cost versus direct accumulation cost
