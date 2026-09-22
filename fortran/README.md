# Fortran backends

`orbital.f90` provides a C-interoperable double-precision magnetic volume kernel for arbitrary orbital positions. The Python package compiles and calls it through `gravmagpy.predict_field` and `gravmagpy.fit_sources`.

The same module exposes `dipole_field` for evaluating fitted equivalent-source moments at arbitrary positions through `predict_equivalent_field`.

The Python model interface (`block`, `polygon`, `equivalent_source_grid`, `sphere_model`, and `run_sphere_model`) generates cards internally and runs `gravmag_sphere_bxyz` or `gravmag_sphere_gauss`. `fit_equivalent_sources` calls `gravmag_sphere_dipole_grid_fit` using normalized observation CSVs. Lower-level `build_fortran`, `run_fortran`, and `run_grid_model` expose compilation, positional arguments, and card execution. The [shell build workflow](../examples/GravMagSphere_README.md) lives in `examples`.

The four executable programs share `gravmag_paths.f90`. Omitted or empty output arguments create files in an `output/` folder beside the input; explicit arguments select a destination. Compile this utility before each executable when building manually. Python and shell builders include it automatically. See [output path conventions](../docs/output_paths.md) for filenames and positional argument examples.
