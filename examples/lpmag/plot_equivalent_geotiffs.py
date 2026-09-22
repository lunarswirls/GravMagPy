"""compare a saved wake-data equivalent-source solution with georeferenced lunar maps"""

from pathlib import Path

from gravmagpy.maps import plot_equivalent_maps
from gravmagpy.utils import artifact_path

# run equivalent_wake.py first to generate this fitted source table
wake_dir = Path("/Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_shadow_state")
source_path = wake_dir / "output" / "equivalent_wake" / "dipoles.csv"
figs_dir = artifact_path(source_path, kind="figs").parent / "equivalent_geotiffs"
map_root = Path("/Users/danywaller/Projects/moon")
bounds = (-108.5, -28.5, -101.5, -21.5)
altitude_km = 30.0
radius_km = 1737.4
shape = (281, 281)

# independent panels use their own units and color limits, like lavapy spectral_plots.py
# replace or extend these dictionaries with lamp reflectance, uv ratio, or other geotiffs
rasters = [
    {"path": map_root / "gradiometry_algorithm/Amitis_test/LROC_WAC_IF_lowres.tif",
     "title": "LROC WAC I/F mosaic", "unit": "I/F", "cmap": "gray", "vmin": 0.03, "vmax": 0.2},
    {"path": map_root / "Magnetic_model_30km_Hood_2022/Hood_2021_modeled_Bmag_30km_65N65S_0.5ppd.tif",
     "title": "Hood magnetic model | 30 km", "unit": "nT", "cmap": "viridis"},
]
if not source_path.is_file():
    raise FileNotFoundError(f"run examples/lpmag/equivalent_wake.py first, or set source_path: {source_path}")
missing = [str(spec["path"]) for spec in rasters if not spec["path"].is_file()]
if missing:
    raise FileNotFoundError(f"set rasters to available geotiff files: {missing}")

image_path = plot_equivalent_maps(
    source_path, rasters, figs_dir / f"equivalent_context_{altitude_km:g}km.png", bounds=bounds,
    altitude_km=altitude_km, radius_km=radius_km, shape=shape, components=("btot",),
    source_layer=1, contours=True, contour_levels=5,
    title="Wake-data equivalent sources and lunar context maps",
)

# a second figure shows all cartesian components alongside the same map products
components_path = plot_equivalent_maps(
    source_path, rasters, figs_dir / f"equivalent_components_{altitude_km:g}km.png", bounds=bounds,
    altitude_km=altitude_km, radius_km=radius_km, shape=shape, components=("bx", "by", "bz", "btot"),
    show_source_panel=False, contours=True, contour_levels=5,
    title="Equivalent-source vector field and lunar context maps",
)
print(f"wrote {image_path}")
print(f"wrote {components_path}")
