# Orientation + Multi-body Example Inputs

Descriptions of gravity and magnetic example input files for GravMag sphere solver (`gravmag_sphere_bxyz`).

All files follow the naming convention: 
`gravmag_sphere_{numbod}body_{ifield}_{nblim}_inc{incdeg}_dec{decdeg}_{opt:desc}`  
where
- `numbod` = total number of sources defined
- `ifield` = one of (mag, grav), describes field controlled by `ifield` flag
- `nblim` = one of (fixedlim, polygon), describes source geometry controlled by `nblim` flag
- `incdeg` = source inclination with respect to Cartesian frame of reference +Z aligned with geographic pole, +X points to sun, +Y eastward
- `decdeg` = source declination with respect to Cartesian frame of reference +Z aligned with geographic pole, +X points to sun, +Y eastward
- `opt:desc` = optional further description of example, i.e., weak magnetization, complex large polygon

## Magnetic 1-body pure orientation examples
- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec90_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc-45_dec120_orient.in`

All use:
- magnetic mode (`ifield=2`)
- fixed body limits (`nblim=1`)
- same geometry/depth, only magnetization orientation differs (`incdeg`, `decdeg`)

## Magnetic multi-body polygon examples
- `gravmag_sphere_3body_mag_polygon_incmix_decmix.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix.in`

Both use:
- repeated body blocks (e.g., `numbod=3`)
- polygon geometry (`nblim=0`)
- mixed orientations and depths

## Magnetic 1-body complex non-rectangular polygon examples
- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in`

This case uses:
- one large irregular polygon with 25 vertices (not rectangular)
- bulk magnetization direction (`incdeg=30`, `decdeg=210`)
- deeper thickness (`depth_top=0.5 km`, `depth_bot=8.0 km`)

## Weak-magnetization examples
- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_weak.in`
- `gravmag_sphere_3body_mag_polygon_incmix_decmix_weak.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix_weak.in`
- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_weak.in`

These examples use the same geometry as their corresponding base cases,
but with reduced magnetization amplitudes to test low-amplitude response.
- magnetization strength (base `M=0.25 A/m`, weak `M=0.05 A/m`)

## Gravity examples
- `gravmag_sphere_1body_grav_fixedlim_incna_decna.in`
- `gravmag_sphere_1body_grav_polygon_incna_decna.in`
- `gravmag_sphere_3body_grav_polygon_incna_decna.in`

These examples use:
- gravity mode (`ifield=1`)
- inclination and declination not applicable ( `incdeg=na`, `decdeg=na`)

## Run examples

```bash
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in output/1body_mag_inc45_dec45.txt
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_3body_mag_polygon_incmix_decmix.in output/3body_mag_incmix_decmix.txt 1
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_5body_mag_polygon_incmix_decmix.in output/5body_mag_incmix_decmix.txt 1
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in output/1body_mag_polygon_complexlarge.txt 1
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_1body_grav_fixedlim_incna_decna.in output/1body_grav_incna_decna.txt
./gravmag_sphere_bxyz 1737.4 examples/gravmag_sphere_3body_grav_polygon_incna_decna.in output/3body_grav_incna_decna.txt 1
```

The optional last argument is `refine_factor`; `1` is faster for multi-body testing.
