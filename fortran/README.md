# Fortran backends

`orbital.f90` provides a C-interoperable double-precision magnetic volume kernel for arbitrary orbital positions. The Python package compiles and calls it through `gravmagpy.predict_field` and `gravmagpy.fit_sources`.

The same module also exposes `dipole_field` for mapping fitted equivalent-source moments at new positions/altitudes through `predict_equivalent_field`. This does not modify the equivalent-grid fitter or any input-card interface.

The Python model interface (`block`, `polygon`, `equivalent_source_grid`, `sphere_model`, and `run_sphere_model`) generates cards internally and runs `gravmag_sphere_bxyz` or `gravmag_sphere_gauss`. `fit_equivalent_sources` calls `gravmag_sphere_dipole_grid_fit` using normalized observation CSVs. Lower-level `build_fortran`, `run_fortran`, and `run_grid_model` remain available. The shell build workflow lives in `examples/gravmag_sphere`.

The four executable programs share `gravmag_paths.f90`. Omitted or empty output arguments now create files in an `output/` folder beside the input; explicit paths and their positional slots are unchanged. Compile this utility before each executable when building manually. Python and shell builders already include it. See [output path conventions](../docs/output_paths.md) for filenames and examples.
