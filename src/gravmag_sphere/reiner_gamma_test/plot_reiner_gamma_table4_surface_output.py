#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path

LOCAL_HERE = Path(__file__).resolve().parent
LOCAL_CACHE = LOCAL_HERE / ".cache"
LOCAL_MPLCONFIG = LOCAL_HERE / ".mplconfig"
LOCAL_CACHE.mkdir(exist_ok=True)
LOCAL_MPLCONFIG.mkdir(exist_ok=True)
os.environ["XDG_CACHE_HOME"] = str(LOCAL_CACHE)
os.environ["MPLCONFIGDIR"] = str(LOCAL_MPLCONFIG)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent

# Hardcoded files and plot settings for the Reiner Gamma Table 4 surface case.
INPUT_IN = HERE / "reiner_gamma_table4_surface_direct.in"
INPUT_BRTP = HERE / "reiner_gamma_table4_surface_direct_brtp.txt"
OUTPUT_PNG = HERE / "reiner_gamma_table4_surface_direct_brtp_2x2.png"
BOX_WIDTH_DEG = 0.0196
THICKNESS_KM = 1.0
ELVO_KM = 0.0
LAT_MIN_DEG = 5.5
LAT_MAX_DEG = 9.0
LON_MIN_E_DEG = 299.5
LON_MAX_E_DEG = 303.5


def load_brtp_sum(path: Path) -> dict[str, np.ndarray]:
    arr = np.loadtxt(path, comments="#")
    if arr.ndim == 1:
        arr = arr[None, :]

    coords = np.column_stack([arr[:, 1], arr[:, 2]])
    uniq, inv = np.unique(coords, axis=0, return_inverse=True)

    br = np.zeros(len(uniq))
    bt = np.zeros(len(uniq))
    bp = np.zeros(len(uniq))
    np.add.at(br, inv, arr[:, 3])
    np.add.at(bt, inv, arr[:, 4])
    np.add.at(bp, inv, arr[:, 5])
    btot = np.sqrt(br * br + bt * bt + bp * bp)

    lon_u = np.unique(uniq[:, 0])
    lat_u = np.unique(uniq[:, 1])

    def reshape_field(values: np.ndarray) -> np.ndarray:
        grid = np.full((len(lat_u), len(lon_u)), np.nan)
        lon_to_ix = {v: i for i, v in enumerate(lon_u)}
        lat_to_iy = {v: i for i, v in enumerate(lat_u)}
        for (lo, la), val in zip(uniq, values):
            grid[lat_to_iy[la], lon_to_ix[lo]] = val
        return grid

    return {
        "lon": lon_u,
        "lat": lat_u,
        "Br": reshape_field(br),
        "Btheta": reshape_field(bt),
        "Bphi": reshape_field(bp),
        "Btot": reshape_field(btot),
    }


def make_field_plot(path: Path, fields: dict[str, np.ndarray], subtitle: str) -> None:
    keys = ["Br", "Btheta", "Bphi", "Btot"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    lon = fields["lon"]
    lat = fields["lat"]

    for ax, key in zip(axes.flat, keys):
        data = fields[key]
        if key == "Btot":
            vmax = float(np.nanmax(data))
            im = ax.imshow(
                data,
                origin="lower",
                aspect="auto",
                extent=[lon.min(), lon.max(), lat.min(), lat.max()],
                cmap="viridis",
                vmin=0.0,
                vmax=vmax,
            )
        else:
            vmax = float(np.nanmax(np.abs(data)))
            im = ax.imshow(
                data,
                origin="lower",
                aspect="auto",
                extent=[lon.min(), lon.max(), lat.min(), lat.max()],
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
            )
        ax.set_title(key)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        fig.colorbar(im, ax=ax, shrink=0.9, label="nT")

    fig.suptitle(f"Reiner Gamma Table 4 direct-solver approximation\n{subtitle}")
    fig.savefig(path, dpi=220)
    plt.close(fig)


if not INPUT_IN.is_file():
    raise FileNotFoundError(f"Missing input file: {INPUT_IN}")

if not INPUT_BRTP.is_file():
    raise FileNotFoundError(
        f"Missing BRTP output: {INPUT_BRTP}. Run run_reiner_gamma_table4_surface_direct.sh first."
    )

fields = load_brtp_sum(INPUT_BRTP)
subtitle = (
    f"55 Table-4 dipoles approximated as {BOX_WIDTH_DEG:.4f} deg square x "
    f"{THICKNESS_KM:.1f} km boxes, altitude {ELVO_KM:.1f} km, "
    f"window {LAT_MIN_DEG:.1f}-{LAT_MAX_DEG:.1f} N / "
    f"{LON_MIN_E_DEG:.1f}-{LON_MAX_E_DEG:.1f} E"
)
make_field_plot(OUTPUT_PNG, fields, subtitle)

print(f"Input IN: {INPUT_IN}")
print(f"Input BRTP: {INPUT_BRTP}")
print(f"Wrote plot: {OUTPUT_PNG}")
