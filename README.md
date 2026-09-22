# GravMagSphere

Python interface to Fortran gravity/magnetic models and subsurface magnetic-source inversion from spacecraft measurements. Python handles CSV ingestion, geometry, bounded minimization, diagnostics, and plots; Fortran evaluates the magnetic forward model.

## Repository layout

```text
src/gravmagpy/              importable Python package
  data.py                  wake-sorted LPMAG CSV ingestion
  geometry.py              spherical block geometry and coordinate conversion
  forward.py               double-precision orbital Fortran interface
  inversion.py             weighted geometry/magnetization fitting and output
  cards.py                 .in cards, grid runs, fitted-model export
  sphere.py                Python block/polygon models, equivalent volume grids, solver calls
  equivalent.py            Python interface to the equivalent dipole-grid fitter
  plotting.py              Bx/By/Bz/Btot grid and orbital-fit plots
  maps.py                  georeferenced GeoTIFF panels and equivalent-source maps
  utils/fortran.py          cached compilation and GravMag Sphere executable calls
  utils/paths.py            input-relative output and figure destinations
fortran/
  orbital.f90              shared-library kernel for arbitrary spacecraft positions
  gravmag_*.f90             grid, spectral, conversion, and dipole-fit programs
examples/                  test cases, notebooks, shell workflows, and saved results
  lpmag/                   real wake-data and synthetic recovery workflows
diagnostics/               benchmark/comparison scripts, external-solver inputs, reports, plots
docs/                      architecture, modeling, mapping, and execution guides
tests/                     numerical and integration checks
data/spice_kernels/         SPICE kernels, separate from package code
```

## Install

Numeric outputs default to `output/` beside the input file, and figures to sibling `figs/`. Explicit destinations also work; see [output path conventions](docs/output_paths.md) for Python, Bash, and Fortran examples.

Use a virtual environment and an installed `gfortran` compiler:

```bash
/Users/danywaller/code/venvs/gravmagpy/bin/python -m pip install -e '.[plot]'
```

Imports do not compile or run models. `build_fortran()` or the first forward call compiles into a temporary-directory cache keyed by source contents, compiler version, flags, and platform. Pass `compiler=` and `build_dir=` to `build_fortran` for other locations. Wheels include the Fortran source and compile locally when used; they do not include machine-specific binaries. See [the architecture notes](docs/architecture.md).

## Fit orbital measurements

```python
from pathlib import Path
from gravmagpy import load_lpmag_csv, fit_sources, save_fit, write_source_input

wake_dir = Path('/Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_shadow_state')
observations = load_lpmag_csv(
    sorted(wake_dir.glob('ma*_nightside.csv')),
    latitude_range=(-27.5, -22.5), longitude_range=(-107.5, -102.5),
    altitude_range_km=(10, 130), wake_only=True, sigma_nt=1.0,
)
sources = [{
    'lat_deg': -25.0, 'lon_deg': -105.0,
    'lat_width_deg': 2.0, 'lon_width_deg': 2.0,
    'depth_top_km': 5.0, 'thickness_km': 5.0,
    'mx_a_m': 1.0, 'my_a_m': 1.0, 'mz_a_m': 1.0,
}]
bounds = [{
    'lat_deg': (-27, -23), 'lon_deg': (-107, -103),
    'depth_top_km': (1, 30),
    'mx_a_m': (-20, 20), 'my_a_m': (-20, 20), 'mz_a_m': (-20, 20),
}]
result = fit_sources(observations, sources, bounds, fit_background=True)
output_dir = wake_dir / 'output' / 'my_fit'
save_fit(result, observations, output_dir)
write_source_input(result['sources'], output_dir / 'sources.in')
```

Each observation retains its actual three-dimensional position and altitude. Omitted parameters stay fixed. Multiple source dictionaries fit multiple blocks together. `fit_background=True` fits one constant vector per file; supply `observations['groups']` for orbit/pass-specific offsets. The result reports convergence, residuals, bound activity, and Jacobian rank; the fitted geometry is conditional on the chosen source count, bounds, magnetization assumptions, and background model.

The runnable [real-data example](examples/lpmag/fit_wake.py) fits location, horizontal dimensions, depth, and magnetization, using the configured wake-data path above. The [synthetic example](examples/lpmag/synthetic_fit.py) generates a CSV with the same SEL column names at 20, 40, and 80 km and recovers a known source:

```bash
PYTHONPATH=src /Users/danywaller/code/venvs/gravmagpy/bin/python examples/lpmag/synthetic_fit.py
PYTHONPATH=src /Users/danywaller/code/venvs/gravmagpy/bin/python examples/lpmag/fit_wake.py
```

See the architecture's [LPMAG reference workflow](docs/architecture.md#lpmag-reference-workflow) for coordinate conventions, quality selection, uncertainties, priors, and practical limits, followed by the path to body-independent, multi-solver equivalent sources.

For geologic context maps, use `gravmagpy.maps.plot_equivalent_maps` to place the fitted equivalent field and source moments alongside GeoTIFF products. The [runnable map example](examples/lpmag/plot_equivalent_geotiffs.py) uses a local WAC mosaic and a 30 km magnetic-model map; [the mapping guide](docs/maps.md) explains CRS handling, spectral-map settings, and prediction at new altitudes. Install `.[maps]` for the optional Rasterio/Matplotlib dependencies.

## Python models for GravMag Sphere

The GravMag Sphere tools are maintained as backends with Python frontends. You can define sources and grids in Python; input cards are generated internally:

```python
import numpy as np
from gravmagpy import equivalent_source_grid, sphere_model, run_sphere_model, write_model_input

sources = equivalent_source_grid(
    latitude_edges_deg=[6, 7, 8], longitude_edges_deg=[-61, -60, -59],
    depth_edges_km=[5, 10], magnetization_a_m=[1, 0, 0],
)
model = sphere_model(sources, latitude_deg=np.linspace(5, 9, 21),
                     longitude_deg=np.linspace(-62, -58, 21), altitude_km=40)
result = run_sphere_model(model, solver="direct", solver_options={"refine_factor": 2})
bx, by, bz, btot = (result[key] for key in ("bx_nt", "by_nt", "bz_nt", "btot_nt"))
write_model_input(model, "output/equivalent_layer.in")
```

Use `block` or `polygon` for individual magnetic or gravity volumes. `fit_equivalent_sources(observations, depth_km=10, spacing_deg=(1, 1))` calls the Fortran dipole-grid fitter, using the same observation dictionary as `load_lpmag_csv`. It estimates point-dipole moments on a fixed grid; `fit_sources` instead optimizes finite-volume geometry. See the [Python modeling guide](docs/sphere.md) for runnable examples, units, regularization, and import/export.

## Fortran input-card models

```python
from gravmagpy import run_grid_model
from gravmagpy.plotting import plot_fields

input_path = 'examples/seafloor_spreading_cases/st_paul_symmetric_stripes.in'
result = run_grid_model(input_path, radius_km=6371.2, solver='direct', options=(1,))
plot_fields(input_path, result['output_path'])
```

`build_fortran` supports `orbital`, `direct`, `spectral`, `gauss_legendre`, `xyz_to_brtp`, and `dipole_grid`. `run_fortran` exposes the positional interfaces. Shell runners work from `examples`; their builds resolve the canonical files directly under `fortran`. Fortran `.in` files are supported directly. See [architecture and development direction](docs/architecture.md) and the [solver guide](examples/GravMagSphere_README.md).

Use `run_sphere_model(model, solver="gauss_legendre")` to evaluate blocks or polygons with the restored Gauss–Legendre volume method in double precision. Node orders and composite subdivisions are configurable; see the [Python modeling guide](docs/sphere.md). The [three-solver diagnostic report](diagnostics/gauss_legendre_comparison/README.md) compares quadrature convergence with direct and pure spectral fields at surface, orbital and far-field altitudes.

## Tests and diagnostics

Solver comparisons and compiler diagnostics are grouped in the top-level [diagnostics folder](diagnostics/README.md). Package regression tests are in `tests/`.

Tests exercise the Fortran kernel against an independent point-dipole limit and the GravMag Sphere surface-charge solver, quadrature convergence, magnetic cancellation, CSV handling, inversion, Python model/card compatibility, gravity and spectral calls, and equivalent-dipole fitting.

## Authors and acknowledgment

Primary author: [Dany Waller](https://danywaller.github.io), [dany.c.waller@gmail.com](mailto:dany.c.waller@gmail.com).

GravMag Sphere credits Dr. Dhahanjay Ravat for the `sphere` program and the University of Kentucky Fall 2018 EES-395 course for its modular formulation. The Fortran sources contain attribution and numerical references.
