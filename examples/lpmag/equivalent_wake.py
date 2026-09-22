"""fit a fixed equivalent dipole layer with the maintained gravmag sphere tool"""

from pathlib import Path

from gravmagpy import load_lpmag_csv, fit_equivalent_sources
from gravmagpy.plotting import plot_fit

wake_dir = Path("/Users/danywaller/Projects/moon/lpmag_l1b_5s_avg/lp_mag_shadow_state")
output_dir = wake_dir / "output" / "equivalent_wake"
figs_dir = wake_dir / "figs" / "equivalent_wake"
paths = sorted(wake_dir.glob("ma*_nightside.csv"))
observations = load_lpmag_csv(
    paths, wake_only=True, sigma_nt=1.0, latitude_range=(-27.5, -22.5),
    longitude_range=(-107.5, -102.5), altitude_range_km=(10.0, 130.0),
)
# this tool fits the supplied vectors without per-file background subtraction
result = fit_equivalent_sources(
    observations, depth_km=10.0, spacing_deg=(1.0, 1.0), layers=1,
    padding_deg=(0.5, 0.5), regularization=0.0, output_dir=output_dir,
)
plot_fit(observations, result, figs_dir / "fit.png")
print(f"fit {len(result['source_table'])} dipoles to {len(observations['xyz_km'])} observations")
print(f"component-averaged rmse: {result['rmse_nt']:.3f} nt")
print(f"wrote {output_dir}")
