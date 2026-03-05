# GravMagSphere Inclination/Declination Physical Audit (2026-03-04)

## Scope
This audit investigates:
- how GravMagSphere interprets inclination/declination,
- whether the implementation is physically coherent,
- how to interpret GravMagSphere vs SHTOOLS comparison residuals.

Primary artifacts audited:
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/gravmag_sphere_bxyz.f90`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/gravmag_sphere_gauss.f90`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/gravmag_sphere_physics.f90`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/gravmag_xyz_to_brtp.f90`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/external_solver_compare.py`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/external_solver_comparison.csv`
- `/Users/danywaller/code/GravMagPy/src/gravmag_sphere/diagnostics/gravmag_vs_shtools_residuals.md`

---

## 1) Coordinate System and Angle Convention Used by GravMagSphere

### 1.1 Global Cartesian frame
Both direct and spectral solvers use the same global planet-fixed Cartesian frame:
- `x = r cos(lat) cos(lon)`
- `y = r cos(lat) sin(lon)`
- `z = r sin(lat)`

Source code:
- `gravmag_sphere_bxyz.f90` lines 348-350, 678-680
- `gravmag_sphere_gauss.f90` lines 493-495, 425-427

Interpretation:
- +X points to lon=0, lat=0.
- +Y points to lon=90, lat=0.
- +Z points to the north pole.

### 1.2 Inclination/declination mapping
Both solvers convert Card-5 `(M_amp, inc, dec)` to Cartesian magnetization as:
- `Mx = M cos(I) cos(D)`
- `My = M cos(I) sin(D)`
- `Mz = M sin(I)`

Source code:
- `gravmag_sphere_bxyz.f90` lines 247-254
- `gravmag_sphere_gauss.f90` lines 322-329

This is mathematically consistent and normalized (unit-direction check below), but it is a **global XYZ angular convention**:
- declination `D` is azimuth in global XY plane from +X toward +Y,
- inclination `I` is elevation from global XY plane toward +Z.

Important caveat:
- This is **not** the usual local geophysical declination/inclination convention (defined in local horizontal N/E plane with local vertical).

---

## 2) Physics Implementation Review

### 2.1 Direct solver magnetic physics (`gravmag_sphere_bxyz`)
Magnetic field is computed with equivalent magnetic surface charge:
- elemental charge: `q = (M · n_out) dS`
- field kernel: `dB = mu0/(4pi) * q * R / |R|^3`

Source code:
- physical kernel: `gravmag_sphere_physics.f90` lines 113-146
- top/bottom/sides source construction: `gravmag_sphere_bxyz.f90` lines 388-477
- side-wall triangle outward normal orientation by interior-point test: lines 755-761
- elemental charge assignment: line 764

Assessment:
- This is a standard and physically valid formulation for uniformly magnetized bodies.

### 2.2 Spectral solver magnetic physics (`gravmag_sphere_gauss`)
Spectral solver first builds volume dipole elements (`m = M dV`) then fits SH coefficients:
- dipole kernel: `B = mu0/(4pi) * (3 (m·R) R / |R|^5 - m / |R|^3)`

Source code:
- source dipoles from `Mx,My,Mz`: `gravmag_sphere_gauss.f90` lines 437-440
- dipole kernel (double precision): lines 1014-1042
- fit samples and SH solve: lines 502-583

Assessment:
- This is physically valid; residuals relative to direct model are dominated by fit/truncation/regularization effects, not by an obvious inc/dec sign/frame bug.

### 2.3 Spherical component conventions
`Br/Btheta/Bphi` conversion is internally consistent:
- solver uses spherical basis with `theta` as colatitude and `Btheta` southward,
- converter uses inverse-consistent transform.

Source code:
- XYZ from `Br/Btheta/Bphi`: `gravmag_sphere_gauss.f90` lines 760-763
- `Br/Btheta/Bphi` from XYZ: `gravmag_xyz_to_brtp.f90` lines 141-143

A direct numerical roundtrip test of these formulas gave max error `1.33e-15`.

---

## 3) Targeted Numerical Validation Checks

### 3.1 Unit-direction check for inc/dec formula
Random Monte Carlo test over `(I,D)` confirmed
`||[cosI cosD, cosI sinD, sinI]|| = 1` to machine precision.
- max unit-vector norm error: `2.22e-16`

### 3.2 Declination invariance at vertical inclination (I=90)
Expected: if `I=90`, `cos(I)=0`, so dec should not matter.

Test (direct solver, same geometry, `dec=0` vs `dec=137`):
- max `|Btot diff| = 7.0e-06 nT`
- RMSE `Btot diff = 3.34e-07 nT`
- relative max diff `1.41e-05 %` of max field

Test (spectral solver, same geometry, `dec=0` vs `dec=137`):
- same chosen auto params in both runs (`lmax=18`, `reg_lambda=0.05`, `reg_power=6`)
- max `|Btot diff| = 2.13e-02 nT`
- RMSE `Btot diff = 1.75e-02 nT`
- relative max diff `0.116 %` of max field

Interpretation:
- Direct path is essentially dec-invariant at `I=90`.
- Spectral path shows tiny leakage from finite-precision + inversion sensitivity, but magnitude is very small.

### 3.3 Orientation sanity at map center (I=0)
Direct-solver center-point fields for fixed-limits orientation test:
- `I=0, D=0`: center field unit direction ~ `[+1, 0, 0]`
- `I=0, D=90`: center field unit direction ~ `[0, -1, 0]`

This confirms declination rotates the effective magnetization azimuth in the global XY frame as implemented.

---

## 4) What SHTOOLS Comparison Does (and Does Not) Test

### 4.1 In the current comparison pipeline, SHTOOLS does not consume inc/dec
`external_solver_compare.py` fits SH coefficients directly to observed component fields (`Br`, `Btheta`, `Bphi`) via LSQ.
- `SHExpandLSQ` and `LSQ_G` are used on component data.
- inc/dec is not read by SHTOOLS in that pipeline.

Source code:
- `solve_shtools_lsq_component`: `external_solver_compare.py` lines 381-399
- case comparison loop: lines 661-789

Implication:
- GravMagSphere vs SHTOOLS residuals mostly diagnose solver approximation/modeling differences, not inc/dec parsing differences.

### 4.2 Residual behavior in latest comparison output
From `diagnostics/external_solver_comparison.csv` (auto-mode GravMagSphere spectral + optimized SHTOOLS/SciPy grids):

- Fixed-limit magnetic orientation cases (non-weak):
  - GravMagSphere mean RMSE(Btot): `2.87`
  - SHTOOLS mean RMSE(Btot): `4.19`
  - SciPy mean RMSE(Btot): `4.26`

- Multi-body polygon cases (non-weak):
  - GravMagSphere mean RMSE(Btot): `4.91`
  - SHTOOLS mean RMSE(Btot): `4.84`
  - SciPy mean RMSE(Btot): `4.88`

- Complex large polygon case:
  - GravMagSphere RMSE(Btot): `318.99`
  - SHTOOLS RMSE(Btot): `278.47`
  - SciPy RMSE(Btot): `279.29`

Interpretation:
- No case pattern suggests an inc/dec sign/axis bug as the dominant residual source.
- Largest errors concentrate in geometrically sharp/complex cases where spectral truncation and ringing/regularization dominate.

---

## 5) Physical Correctness Judgment

### Verdict
The inclination/declination implementation is **physically coherent and internally consistent** with the model's global XYZ coordinate system.

What is correct:
- same `inc/dec -> Mx,My,Mz` mapping in direct and spectral paths,
- physically valid magnetic kernels in both paths,
- consistent Cartesian/spherical transforms,
- no evidence of a direct sign inversion or axis swap bug.

What can be misinterpreted:
- The angle convention is global XYZ, not local geophysical dec/inc.
- If users assume local N/E/U (or N/E/Down) conventions, they can interpret inputs incorrectly even though solver math is consistent.

---

## 6) Recommended Improvements

1. Explicitly document angle convention in README and input-card docs.
- State that Card-5 declination is azimuth from global +X toward +Y; inclination is from global XY toward +Z.

2. Add optional local geophysical mode.
- Add a switch to interpret dec/inc in local frame at body center (N/E/U or N/E/Down), then rotate to global XYZ.

3. Add a small numerical guard at `|I| ~ 90`.
- If `abs(cosI) < eps`, set horizontal components exactly to zero to suppress tiny dec leakage in spectral fitting.

4. Add regression tests.
- Vertical-inclination dec-invariance test (`I=90`, vary `D`).
- Known-axis tests (`I=0,D=0`; `I=0,D=90`; `I=90,D=any`).
- Cartesian/spherical roundtrip consistency tests.


## 9) Implications for Cross-Model Residual Interpretation

Given the audit findings, model differences should be interpreted as a combination of:

- physically different inverse assumptions (equivalent-source versus pure SH LSQ formulations),
- coordinate-convention differences (global XYZ angle interpretation versus local geophysical conventions in other tools),
- regularization and bandwidth choices that redistribute anomaly energy.

Practical consequence:

- residuals near sharp edges mainly indicate basis/regularization limits,
- broad regional offsets indicate different low-degree parameterization or continuation behavior,
- multi-body residuals indicate coefficient competition and source-partition choices rather than a single coding bug.
