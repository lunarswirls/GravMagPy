# Reiner Gamma Test
## Hemingway & Garrick-Bethell Table 4 definitions
- See Hemingway and Garrick-Bethell (2012), "Magnetic field direction and lunar swirl morphology: Insights from Airy and Reiner Gamma" :)
- Table columns: `['Latitude', 'Longitude', 'Depth', 'Magnetic Moment', 'Inclination', 'Declination']`
- Units: `['deg N', 'deg E', 'km', 'A-m^2', 'deg', 'deg']`
- Number of dipoles: `55`
- Depth from Table 4: `5.0 km`
- Equal moment per dipole: `1.820e+11 Am2`
- Local dipole orientation: `I=2.0 deg`, `D=-8.0 deg`

## GravMagPy definitions
- Global GravMagPy orientation: `inc=77.639 deg`, `dec=161.754 deg`
  - Obtained by converting the Table 4 local direction to global direction at the source-set centroid (`lat=7.439 deg`, `lon=-58.765 deg`).
- Finite-body approximation: `0.0196 deg x 0.0196 deg x 1.0 km`
  - This is necessary because Table 4 defines dipole-array source centers, while `gravmag_sphere_bxyz` models finite uniformly magnetized bodies.
- Constant magnetization strength derived from the Table-4 dipole moment using the centroid box volume, then applied to every source body.
- Uniform magnetization in Card 5: `5.248018e+02 A/m`
  - Because the finite boxes sit at slightly different latitudes, their physical volumes differ slightly and therefore their implied finite-body moments vary slightly even though all Card-5 magnetization strengths are identical.
- Evaluation window: `5.5 to 9.0 deg N`, `299.5 to 303.5 deg E` (`0 to 360` convention)
  - Wrapped solver window: `-60.5 to -56.5 deg` (`-180 to 180` convention)
- Observation grid: `dlat=dlon=0.020`, `nlat=176`, `nlon=201`
