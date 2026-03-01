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
cd /Users/tycho/code/GravMagPy/src/gravmag_sphere
./build_gravmag_tools.sh
```

Build outputs:

- `./gravmag_sphere_brtp`
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
  [source_nlat] [source_nlon] [source_nr]
```

Example:

```bash
./gravmag_sphere_gauss 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in \
  output/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base_xyz.txt \
  24 2 72 144 0.2 4.0 0 0 0
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
  [source_nlat] [source_nlon] [source_nr]
```

Examples:

```bash
./run_input_to_xyz.sh direct 1737.4 \
  examples/gravmag_sphere_1body_grav_fixedlim_incna_decna.in

./run_input_to_xyz.sh spectral 1737.4 \
  examples/gravmag_sphere_1body_mag_fixedlim_inc90_dec0_base.in
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
./run_gravmag_sphere_gauss.sh <R_sphere_km> <input.in> [xyz_output] [lmax] [refine] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] [source_nlat] [source_nlon] [source_nr]
```

## Visualization scripts

### 2D maps for all examples (Python)

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
