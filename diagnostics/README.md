# GravMag Sphere diagnostics

This folder groups compiler benchmarks, direct/spectral comparisons, equivalent-source comparisons, their formatted inputs, and saved reports/figures. I AM ACTIVELY DEBUGGING THE CODE!!

## Contents

- `external_solver_inputs/`: sixteen formatted comparison cases, observations, metadata, and saved solver outputs
- `diagnostic_fortran.py`: compiler benchmarking and spectral-degree ringing sweeps
- `external_solver_compare.py`: external input generation and SciPy/SHTOOLS/Harmonica comparisons
- `hybrid_vertex_sweep.py`: polygon-complexity sweep
- `shtools_operating_map.py`: accuracy/runtime operating maps
- `complexlarge_direct_vs_spectral_analysis.py`: boundary-visibility analysis
- `plot_fixed_polygon_compare.py`: fixed-limit/polygon field comparison plot
- `run_comparison_tests.sh`: comparison-suite runner
- Markdown reports and residual/side-by-side plot directories: existing diagnostic results

Relative paths for benchmark, comparison, and sweep options are resolved from the repository root. The fixed/polygon plotting script uses caller-relative input paths and defaults its output to this folder. The suite runner can be invoked from any directory. 

## Environment and usage

The core diagnostics need NumPy, SciPy, and Matplotlib; SHTOOLS/Harmonica suites additionally need `pyshtools`/`harmonica`. The runner skips missing optional backends unless `--strict-deps` is supplied. Full benchmark/sweep runs rebuild executables and can take substantial time; `--help` is non-mutating.

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/external_solver_compare.py --help
bash diagnostics/run_comparison_tests.sh --help
```

## Workflows

### 1) Compiler runtime benchmark + hotspot report

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/diagnostic_fortran.py benchmark \
  --compilers gfortran,ifx,ifort,flang-new,nvfortran \
  --repeats 3
```

Outputs:

- `diagnostics/fortran_benchmark_report.json`
- `diagnostics/fortran_benchmark_report.md`

The benchmark report includes:

- wall time per compiler/case
- slowest solver stage (from internal Fortran timers)
- largest memory bucket estimate

### 2) Spherical harmonic `lmax` ringing sweep

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python diagnostics/diagnostic_fortran.py lmax-sweep \
  --compiler gfortran \
  --lmin 8 --lmax 96 --lstep 8
```

Outputs:

- `diagnostics/lmax_sweep.csv`
- `diagnostics/lmax_sweep.md`

The sweep compares spectral results to a direct-solver baseline and reports
edge-focused ringing metrics (`edge_amp`, `edge_overshoot`) per `lmax`.

### 3) Direct per-run diagnostics from Fortran executables

Set `GRAVMAG_DIAGNOSTICS=1` to print stage timings and memory estimates:

```bash
GRAVMAG_DIAGNOSTICS=1 examples/gravmag_sphere/gravmag_sphere_gauss 1737.4 examples/gravmag_sphere/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in diagnostics/diag_gauss_xyz.txt 24 2 72 144 0.2 4.0 0 0 0
GRAVMAG_DIAGNOSTICS=1 examples/gravmag_sphere/gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere/lunar_examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in diagnostics/diag_direct_xyz.txt 2 0 0 0
```

Diagnostics lines are emitted as `DIAG|<solver>|<category>|<key>=<value>`.

### 4) Run all solver comparison suites

`run_comparison_tests.sh` orchestrates all comparison workflows:

- external solver comparison (Fortran spectral vs SciPy SH LSQ vs SHTOOLS LSQ)
- GravMagSphere-vs-SHTOOLS residual suite
- GravMagSphere-vs-Harmonica XYZ residual suite
- complexlarge direct-vs-spectral boundary analysis
- hybrid-vs-spectral vertex sweep

```bash
# Run all comparison suites
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python

# Run selected suites only
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python \
  --tests external,shtools,harmonica

# Limit to specific lunar_examples
bash diagnostics/run_comparison_tests.sh --python /Users/danywaller/code/venvs/gravmagpy/bin/python \
  --tests external,harmonica \
  --lunar_examples gravmag_sphere_1body_mag_polygon_inc90_dec0_base,gravmag_sphere_3body_mag_polygon_incmix_decmix
```

Notes:
- Set `--no-clean` to keep intermediate txt/csv artifacts
- Set `--strict-deps` to fail immediately if Python dependencies are missing
