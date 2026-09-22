"""fit spherical subsurface blocks to wake-sorted lpmag observations"""

from pathlib import Path

from gravmagpy import fit_sources, load_lpmag_csv, save_fit, write_source_input
from gravmagpy.plotting import plot_fit

# edit these settings for the anomaly and source hypothesis being tested
# directory containing wake-selected magnetic csv files
wake_dir = Path("/Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_shadow_state")
output_dir = wake_dir / "output" / "wake_fit"
figs_dir = wake_dir / "figs" / "wake_fit"
pattern = "ma*_nightside.csv"
latitude_range = (-27.5, -22.5)
longitude_range = (-107.5, -102.5)
altitude_range_km = (10.0, 130.0)
radius_km = 1737.4

# magnetization is in sel cartesian a/m; depths are below the lunar surface
sources = [{
    "lat_deg": -25.0, "lon_deg": -105.0, "lat_width_deg": 2.0, "lon_width_deg": 2.0,
    "depth_top_km": 5.0, "thickness_km": 5.0, "mx_a_m": 1.0, "my_a_m": 1.0, "mz_a_m": 1.0,
}]
bounds = [{
    "lat_deg": (-27.0, -23.0), "lon_deg": (-107.0, -103.0),
    "lat_width_deg": (0.5, 4.0), "lon_width_deg": (0.5, 4.0), "depth_top_km": (1.0, 30.0),
    "mx_a_m": (-20.0, 20.0), "my_a_m": (-20.0, 20.0), "mz_a_m": (-20.0, 20.0),
}]
# thickness stays fixed to reduce the thickness/magnetization tradeoff
paths = sorted(wake_dir.glob(pattern))
if not paths:
    raise FileNotFoundError(f"no {pattern} files in {wake_dir}")
observations = load_lpmag_csv(
    paths, radius_km=radius_km, wake_only=True, sigma_nt=1.0,
    latitude_range=latitude_range, longitude_range=longitude_range, altitude_range_km=altitude_range_km,
)
print(f"loaded {len(observations['table'])} wake observations from {len(paths)} files")
result = fit_sources(observations, sources, bounds, fit_background=True, quadrature_order=8, max_nfev=150)
save_fit(result, observations, output_dir)
plot_fit(observations, result, figs_dir / "fit.png")
write_source_input(result["sources"], output_dir / "fitted_sources.in", radius_km=radius_km)
print(f"success={result['success']}: {result['message']}")
print(f"component-averaged rmse: {result['initial_rmse_nt']:.3f} -> {result['rmse_nt']:.3f} nt")
print(f"background-only baseline rmse: {result['background_only_rmse_nt']:.3f} nt")
active = [name for name, flag in zip(result["parameter_names"], result["active_bounds"]) if flag]
print(f"parameters at bounds: {active}")
print(f"wrote {output_dir}")
