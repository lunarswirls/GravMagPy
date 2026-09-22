# SHTOOLS and GravMag Sphere spectral comparison

## Scope

The [external-solver report](external_solver_comparison.md) compares Fortran spectral predictions, SciPy spherical-harmonic least squares, and SHTOOLS least squares against direct-solver fields. The [pairwise residual report](gravmag_vs_shtools_residuals.md) compares GravMag Sphere spectral and SHTOOLS predictions directly.

Each report records its solver settings and per-case residuals. Magnetic quantities use nT and gravity quantities use mGal; do not combine their dimensional errors into one aggregate score.

## Model and fitting differences

- GravMag Sphere constructs physical source elements from the input geometry, samples their fields, and estimates spherical-harmonic coefficients with joint-component fitting and degree-dependent regularization
- SciPy and SHTOOLS fit spherical-harmonic expansions to sampled field components independently in these diagnostic adapters
- The Fortran spectral workflow can also apply local edge correction and direct/spectral hybrid evaluation
- Search grids, automatic parameter selection, damping, and component coupling differ between workflows, so equal nominal harmonic degree does not define an equivalent inverse problem

Independent component fits are useful approximation baselines, but are not automatically a physically coupled magnetic equivalent-source solution. The package's fixed-grid dipole inversion is a separate workflow described in [Python modeling](../docs/sphere.md#fit-an-equivalent-point-dipole-model).

## Reading the comparisons

Residuals against the direct solver measure disagreement with that discretized baseline, not error against an exact physical solution. Baseline refinement and independent kernel checks are needed before attributing all disagreement to a spectral method.

Finite harmonic bandwidth limits the representation of sharp spatial variations near source boundaries. Regularization changes the balance of low- and high-degree structure, while local correction and hybrid settings alter near-boundary predictions. Compare component errors as well as total magnitude, because a magnitude match can hide directional disagreement.

Multi-body fields add linearly for prescribed magnetization. Differences in joint fitting, regularization, and source sampling can change approximation errors without implying nonlinear magnetic superposition. Total-field magnitude and source-geometry optimization are nonlinear operations.

Use the [residual methodology](external_solver_residuals_method.md) for grid alignment, residual signs, metrics, and commands. The [architecture roadmap](../docs/architecture.md#path-forward-body-independent-multi-solver-equivalent-sources) describes the operator and validation contracts required for interchangeable production backends.
