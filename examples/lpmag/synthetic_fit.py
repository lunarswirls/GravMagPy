"""create a reproducible lpmag-format csv and recover its source geometry"""

from pathlib import Path

import numpy as np
import pandas as pd

from gravmagpy import fit_sources, load_lpmag_csv, predict_field, save_fit, spherical_to_cartesian, write_source_input
from gravmagpy.plotting import plot_fit

output_dir = Path(__file__).resolve().parent / "output" / "synthetic"
figs_dir = Path(__file__).resolve().parent / "figs" / "synthetic"
output_dir.mkdir(parents=True, exist_ok=True)
radius_km = 1737.4
truth = [{
    "lat_deg": 7.5, "lon_deg": -59.0, "lat_width_deg": 1.2, "lon_width_deg": 1.5,
    "depth_top_km": 8.0, "thickness_km": 6.0, "mx_a_m": 1.5, "my_a_m": -0.7, "mz_a_m": 0.9,
}]
latitude, longitude = np.meshgrid(np.linspace(5.5, 9.5, 13), np.linspace(-61, -57, 11), indexing="ij")
latitude = np.tile(latitude.ravel(), 3)
longitude = np.tile(longitude.ravel(), 3)
altitude = np.repeat([20.0, 40.0, 80.0], 143)
xyz = spherical_to_cartesian(latitude, longitude, radius_km + altitude)
# all measurements are synthetic, noise-free, and labeled explicitly
field = predict_field(xyz, truth, radius_km=radius_km, quadrature_order=10)
table = pd.DataFrame({"utc": pd.date_range("1999-01-01", periods=len(xyz), freq="5s").astype(str)})
for index, key in enumerate(("X_SEL", "Y_SEL", "Z_SEL")):
    table[key] = xyz[:, index]
for index, key in enumerate(("Bx_SEL", "By_SEL", "Bz_SEL")):
    table[key] = field[:, index]
table["Brms"] = 0.0
table["altitude_km"] = altitude
csv_path = output_dir / "synthetic_nightside_wake.csv"
table.to_csv(csv_path, index=False)
observations = load_lpmag_csv(csv_path, radius_km=radius_km)
initial = [dict(truth[0], lat_deg=7.7, lon_deg=-58.8, lat_width_deg=1.0, lon_width_deg=1.8, depth_top_km=12.0)]
bounds = [{
    "lat_deg": (6.0, 9.0), "lon_deg": (-60.5, -57.5), "lat_width_deg": (0.5, 2.0),
    "lon_width_deg": (0.5, 2.5), "depth_top_km": (1.0, 25.0),
}]
result = fit_sources(observations, initial, bounds, quadrature_order=10)
save_fit(result, observations, output_dir)
plot_fit(observations, result, figs_dir / "fit.png")
write_source_input(result["sources"], output_dir / "fitted_sources.in", radius_km=radius_km)
pd.DataFrame(truth).to_csv(output_dir / "truth.csv", index=False)
print(f"success={result['success']}, rmse={result['rmse_nt']:.6g} nt, evaluations={result['nfev']}")
print(f"wrote {output_dir}")
