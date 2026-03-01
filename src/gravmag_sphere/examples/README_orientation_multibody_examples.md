# Orientation + Multi-body Example Inputs

New magnetic examples for the volume solver (`gravmag_sphere_brtp`):

## 1-body orientation series (fixed limits)
- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec0_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc0_dec90_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in`
- `gravmag_sphere_1body_mag_fixedlim_inc-45_dec120_orient.in`

All use:
- `ifield=2` (magnetic mode)
- same geometry/depth, only magnetization orientation differs (`inc`, `dec`)

## Multi-body polygon examples
- `gravmag_sphere_3body_mag_polygon_incmix_decmix.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix.in`
- `gravmag_sphere_3body_mag_polygon_incmix_decmix_weak.in`
- `gravmag_sphere_5body_mag_polygon_incmix_decmix_weak.in`

Both use:
- repeated body blocks
- polygon geometry (`nblim=0`)
- mixed amplitudes, orientations, and depths

## Large complex non-rectangular polygon (1 body)
- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in`
- `gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge_weak.in`

This case uses:
- one large irregular polygon with 25 vertices (not rectangular)
- magnetic mode (`ifield=2`) with `inc=30`, `dec=210`
- deeper thickness (`depth_top=0.5 km`, `depth_bot=8.0 km`)

## Weak-magnetization baseline cases
- `gravmag_sphere_1body_mag_fixedlim_inc90_dec0_weak.in`

These weak examples use the same geometry as their corresponding strong cases,
but with reduced magnetization amplitudes to test low-amplitude response.

## Run examples

```bash
./gravmag_sphere_brtp 1737.4 examples/gravmag_sphere_1body_mag_fixedlim_inc45_dec45_orient.in output/1body_mag_inc45_dec45.txt
./gravmag_sphere_brtp 1737.4 examples/gravmag_sphere_3body_mag_polygon_incmix_decmix.in output/3body_mag_incmix_decmix.txt 1
./gravmag_sphere_brtp 1737.4 examples/gravmag_sphere_5body_mag_polygon_incmix_decmix.in output/5body_mag_incmix_decmix.txt 1
./gravmag_sphere_brtp 1737.4 examples/gravmag_sphere_1body_mag_polygon_inc30_dec210_complexlarge.in output/1body_mag_polygon_complexlarge.txt 1
```

The optional last argument is `refine_factor`; `1` is faster for multi-body testing.
