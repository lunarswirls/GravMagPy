# GravMag Sphere

Fortran + Python workflow for spherical gravity/magnetic forward modeling.

- Fortran solvers produce cartesian field components (`Fx Fy Fz`) on a lon/lat grid
- A Fortran converter rotates cartesian components (`Bx By Bz`) to spherical components (`Br Btheta Bphi`)
- Python scripts create 2D PNG maps and 3D Plotly HTML views

## Directory layout

- `examples/`: input card files (`*.in`)
- `output/`: numeric outputs (`*_xyz.txt`, `*_brtp.txt`)
- `figs/`: figures (`*.png`, `*.html`)

## Build

```bash
cd /Users/danywaller/code/GravMagPy/src/gravmag_sphere
./build_gravmag_tools.sh
```

Build outputs:

- `./gravmag_sphere_bxyz`
- `./gravmag_sphere_gauss`
- `./gravmag_xyz_to_brtp`

## Compiled executables

### 1) Direct volume solver (Fortran)

`gravmag_sphere_bxyz` computes XYZ fields directly from discretized source volume/surfaces

```bash
./gravmag_sphere_bxyz <R_sphere_km> <input.in> <output_xyz.txt> [refine_factor] [source_nlat] [source_nlon] [source_nr]
```

Example:

```bash
./gravmag_sphere_bxyz 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt \
  2 0 0 0
```

### 2) Spectral solver (Fortran)

`gravmag_sphere_gauss` fits spherical harmonic coefficients and evaluates XYZ fields

```bash
./gravmag_sphere_gauss <R_sphere_km> <input.in> <output_xyz.txt> \
  [lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] \
  [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg]
```

Defaults:
- `auto_mode=1` (pilot sweep auto-selects `lmax/reg_lambda/reg_power`)
- `joint_strength=1.0` (joint Br/Btheta/Bphi fit)
- `edge_correction=1` (local edge RBF correction)
- `hybrid_mode=1` (auto: enable hybrid for complex or multi-body inputs)
- `hybrid_band_deg=1.5`
- `complex_vertex_threshold=12`
- `hybrid_transition_deg=0.75`

Example:

```bash
./gravmag_sphere_gauss 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt \
  24 2 72 144 0.2 4.0 0 0 0 1 1.0 1 1 1.5 12 0.75
```

### Spectral Flags and Use Cases

`gravmag_sphere_gauss` / `run_input_to_xyz.sh` / `run_gravmag_sphere_gauss.sh` now expose:

- `auto_mode` (`0|1`):
  - `1` auto-selects `(lmax, reg_lambda, reg_power)` from a pilot sweep.
  - `0` uses user-provided seed values directly.
- `joint_strength` (`>=0`):
  - weight for joint vector fit (`Br/Btheta/Bphi`) in coefficient solve.
  - `0` reduces to radial-only fit.
- `edge_correction` (`0|1`) and `--no-local-correction`:
  - enable/disable local edge RBF correction.
- `hybrid_mode` (`0|1|2`):
  - `0`: spectral-only evaluation.
  - `1`: auto hybrid (enabled when input is complex or multi-body).
  - `2`: force hybrid for all inputs.
- `hybrid_band_deg`:
  - distance from polygon edge where direct kernel result dominates.
- `hybrid_transition_deg`:
  - blend width between direct and spectral predictions.
- `complex_vertex_threshold`:
  - polygon vertex threshold used by complexity classifier in hybrid auto mode.

Typical use cases:

```bash
# 1) Simple body, spectral-only (disable hybrid + edge correction)
./run_input_to_xyz.sh spectral 1737.4 examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  output/simple_spectral_only_xyz.txt \
  --auto-mode 0 --lmax 24 --reg-lambda 0.2 --reg-power 4 \
  --edge-correction 0 --hybrid-mode 0

# 2) Complex polygon, default robust mode (auto + local correction + hybrid auto)
./run_input_to_xyz.sh spectral 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  output/complex_auto_hybrid_xyz.txt

# 3) Force hybrid and widen direct edge band
./run_input_to_xyz.sh spectral 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  output/complex_force_hybrid_xyz.txt \
  --hybrid-mode 2 --hybrid-band-deg 2.0 --hybrid-transition-deg 1.0

# 4) Multi-body collective fit (automatic when bodies share Card-2 grid)
./run_input_to_xyz.sh spectral 1737.4 examples/gravmag_sphere_3body_mag_polygon_incmix_decmix.in \
  output/multibody_collective_xyz.txt --auto-mode 1
```

### 3) XYZ -> spherical converter (Fortran)

`gravmag_xyz_to_brtp` converts XYZ components to `Br Btheta Bphi`

```bash
./gravmag_xyz_to_brtp <input_xyz.txt> <output_brtp.txt>
```

Example:

```bash
./gravmag_xyz_to_brtp \
  output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt \
  output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_brtp.txt
```

## Batch shell scripts

### Single-step wrappers

1. `run_input_to_xyz.sh` (solver step)

```bash
./run_input_to_xyz.sh <solver> <R_sphere_km> <input.in> [xyz_output] \
  [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] \
  [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg] \
  [--lmax N] [--reg-lambda X] [--reg-power P] [--auto-mode 0|1] [--joint-strength X] \
  [--edge-correction 0|1] [--hybrid-mode 0|1|2] [--hybrid-band-deg X] \
  [--complex-vertex-threshold N] [--hybrid-transition-deg X] [--no-local-correction]
```

Examples:

```bash
./run_input_to_xyz.sh direct 1737.4 \
  examples/gravmag_sphere_1body_grav_fixedlim_incna_decna.in

./run_input_to_xyz.sh spectral 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in

./run_input_to_xyz.sh spectral 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  --lmax 30 --reg-lambda 0.1 --reg-power 4 --hybrid-mode 0 --no-local-correction
```

2. `run_xyz_to_brtp.sh` (conversion step)

```bash
./run_xyz_to_brtp.sh <input_xyz.txt> [output_brtp.txt]
```

### End-to-end wrappers

1. `run_gravmag_end_to_end.sh` (single case, direct solver)

```bash
./run_gravmag_end_to_end.sh <R_sphere_km> <input.in> [xyz_out] [brtp_out] \
  [refine_factor] [source_nlat] [source_nlon] [source_nr]
```

2. `run_all_examples_end_to_end.sh` (all examples)

```bash
./run_all_examples_end_to_end.sh [solver] [R_sphere_km] [examples_dir] [output_dir] \
  [refine_factor] [lmax] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] \
  [source_nlat] [source_nlon] [source_nr]
```

Example:

```bash
./run_all_examples_end_to_end.sh spectral 1737.4 examples output 2 24 72 144 0.2 4.0 0 0 0
```

### Solver-specific wrappers

```bash
./run_gravmag_sphere_f90.sh <R_sphere_km> <input.in> [xyz_output] [refine] [source_nlat] [source_nlon] [source_nr]
./run_gravmag_sphere_gauss.sh <R_sphere_km> <input.in> [xyz_output] [lmax] [refine] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] [hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg] [--lmax N] [--reg-lambda X] [--reg-power P] [--auto-mode 0|1] [--joint-strength X] [--edge-correction 0|1] [--hybrid-mode 0|1|2] [--hybrid-band-deg X] [--complex-vertex-threshold N] [--hybrid-transition-deg X] [--no-local-correction]
```

## Visualization scripts

### 2D maps for all examples

`run_all_examples_brtp.py` can run solver + converter + plotting in one Python workflow

```bash
.venv/bin/python run_all_examples_brtp.py \
  --solver spectral \
  --rsphere-km 1737.4 \
  --output-dir output \
  --figs-dir figs
```

Outputs:

- numeric tables in `output/`
- PNG figures in `figs/` (`*_brtp_2x2.png`)

### 3D Plotly for one source case

```bash
.venv/bin/python plot_single_case_3d_plotly.py \
  examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in \
  output/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_brtp.txt \
  gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_3d.html \
  --sparse-step 8 --max-vectors 450 --vector-scale 0.06
```

If the html filename has no directory component, it is written to `figs/`

### Fixed-limit vs polygon comparison plot

```bash
.venv/bin/python plot_fixed_polygon_compare.py \
  --fixed output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_brtp.txt \
  --polygon output/gravmag_sphere_1body_mag_polygon_inc90_dec0_base_brtp.txt \
  --out figs/fixed_vs_polygon_components.png
```

### Helper API for custom plotting

Use `gravmag_plot_helpers.py` from scripts or notebook:

```python
from pathlib import Path
from gravmag_plot_helpers import plot_example

fig = plot_example(
    "gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base",
    root_dir=Path("."),
    examples_dir="examples",
    output_dir="output",
)
fig.savefig("figs/custom_example.png", dpi=250)
```

## Notebook

Open:

- `gravmag_examples_visualization.ipynb`

The notebook is organized as:

1. setup/import cell
2. one markdown explanation + one plotting cell per example
3. each plot cell renders Br/Btheta/Bphi/Btot for that case

Run with Jupyter from this directory:

```bash
.venv/bin/python -m jupyter notebook gravmag_examples_visualization.ipynb
```

## Typical workflow

```bash
# 1) build
./build_gravmag_tools.sh

# 2) run all examples (numeric output)
./run_all_examples_end_to_end.sh spectral

# 3) generate figures
.venv/bin/python run_all_examples_brtp.py --solver spectral --output-dir output --figs-dir figs
```

## Diagnostics

### 1) Compiler runtime benchmark + hotspot report

```bash
python3 diagnostic_fortran.py benchmark \
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
python3 diagnostic_fortran.py lmax-sweep \
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
GRAVMAG_DIAGNOSTICS=1 ./gravmag_sphere_gauss 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in output/diag_gauss_xyz.txt 24 2 72 144 0.2 4.0 0 0 0
GRAVMAG_DIAGNOSTICS=1 ./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc90_dec0_base.in output/diag_direct_xyz.txt 2 0 0 0
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
./run_comparison_tests.sh

# Run selected suites only
./run_comparison_tests.sh \
  --tests external,shtools,harmonica

# Limit to specific examples
./run_comparison_tests.sh \
  --tests external,harmonica \
  --examples gravmag_sphere_1body_mag_polygon_inc90_dec0_base,gravmag_sphere_3body_mag_polygon_incmix_decmix
```

Notes:
- Set `--no-clean` to keep intermediate txt/csv artifacts
- Set `--strict-deps` to fail immediately if Python dependencies are missing
