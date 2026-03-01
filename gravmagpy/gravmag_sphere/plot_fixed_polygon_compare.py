#!/usr/bin/env python3
"""Compare fixed-limit and polygon outputs in a side-by-side 2x4 panel."""

import argparse
import numpy as np
import matplotlib.pyplot as plt


def _parse_header(path):
    """Return first comment-header tokens if present."""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                return s[1:].strip().split()
            break
    return None


def _label_and_unit(name):
    """Split a column label like 'Br_nT' into ('Br','nT')."""
    if "_" not in name:
        return name, ""
    stem, unit = name.rsplit("_", 1)
    return stem, unit


def load_fields(path):
    """Load one solver output file and return named component arrays."""
    header = _parse_header(path)
    a = np.loadtxt(path, comments="#")
    if a.ndim != 2 or a.shape[1] < 7:
        raise ValueError(f"Expected >=7 columns in {path}")

    if header is not None and len(header) >= 7:
        c1, c2, c3, c4 = header[3], header[4], header[5], header[6]
    else:
        c1, c2, c3, c4 = "Bx_nT", "By_nT", "Bz_nT", "Btot_nT"

    fields = {
        "lon": a[:, 1],
        "lat": a[:, 2],
        c1: a[:, 3],
        c2: a[:, 4],
        c3: a[:, 5],
        c4: a[:, 6],
    }
    return fields, [c1, c2, c3, c4]


def to_grid(lon, lat, values):
    """Convert lon/lat/value rows into a dense regular grid."""
    lon_r = np.round(lon, 6)
    lat_r = np.round(lat, 6)
    ulon = np.unique(lon_r)
    ulat = np.unique(lat_r)
    ulon.sort()
    ulat.sort()

    ilon = {x: i for i, x in enumerate(ulon)}
    ilat = {y: i for i, y in enumerate(ulat)}

    z = np.full((ulat.size, ulon.size), np.nan)
    for x, y, vv in zip(lon_r, lat_r, values):
        z[ilat[y], ilon[x]] = vv

    if np.isnan(z).any():
        raise RuntimeError("Grid contains missing cells; inputs must be regular lon/lat grid outputs.")

    return ulon, ulat, z


def plot_compare(fixed, poly, cols, out_png, dpi=250):
    """Plot fixed-limit (row 1) and polygon (row 2) for each component column."""
    label0, unit0 = _label_and_unit(cols[0])
    label1, unit1 = _label_and_unit(cols[1])
    label2, unit2 = _label_and_unit(cols[2])
    label3, unit3 = _label_and_unit(cols[3])

    labels = [label0, label1, label2, label3]
    units = [unit0, unit1, unit2, unit3]
    unit = units[0] if units[0] else units[3]

    fig, axes = plt.subplots(2, 4, figsize=(16, 7), constrained_layout=True)

    for col, key in enumerate(cols):
        for row, data in enumerate((fixed, poly)):
            ax = axes[row, col]
            ulon, ulat, z = to_grid(data["lon"], data["lat"], data[key])

            if col == 3:
                vmin, vmax = 0.0, float(np.nanmax(z))
                cmap = "viridis"
            else:
                m = float(np.nanmax(np.abs(z)))
                vmin, vmax = -m, m
                cmap = "RdBu_r"

            im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("Longitude [deg]")
            if col == 0:
                row_label = "Fixed-limit" if row == 0 else "Polygon"
                ax.set_ylabel(f"Latitude [deg]\\n{row_label}")
            ax.set_title(labels[col])

            cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.8)
            cb.set_label(unit)

    fig.suptitle("Fortran Volume Solver: Fixed-limit vs Polygon")
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)


def main():
    """CLI entrypoint."""
    p = argparse.ArgumentParser(description="Plot fixed-limit vs polygon outputs")
    p.add_argument("--fixed", required=True, help="Fortran output txt for fixed-limit case")
    p.add_argument("--polygon", required=True, help="Fortran output txt for polygon case")
    p.add_argument("--out", default="figs/fixed_vs_polygon_components.png")
    p.add_argument("--dpi", type=int, default=250)
    args = p.parse_args()

    fixed, cols_fixed = load_fields(args.fixed)
    poly, cols_poly = load_fields(args.polygon)

    if cols_fixed != cols_poly:
        raise ValueError(f"Column mismatch between files: {cols_fixed} vs {cols_poly}")

    plot_compare(fixed, poly, cols_fixed, args.out, dpi=args.dpi)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
