#!/usr/bin/env python3
"""
Reads parameters from gravmag_sphere_brtp input file and plots results.txt

Author: Dany Waller
Submitted: Nov 08 2018
Title: EES-395 Final project

Modifications:
2026-02-28
- Updated to support new multi-column results format that outputs:
       lon lat Br Btheta Bphi Btot
    and plots a 2x2 panel: Br, Btheta, Bphi, |B|
2025-01-19
- Parses which field is computed from Card 5: ifield (total / Br / Bphi / Btheta / gravity)
2020-09-28
- Removed PyGMT, converted to Matplotlib
"""
import os
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# edit me! :)
# ============================================================
# input_in_file = "gravmag_sphere_test.in"    # Fortran input file
# results_brtp   = "results.txt"              # Fortran output

input_in_file = "/Users/danywaller/Projects/moon/equivalent_source_modeling/reiner_gamma_3body.in"
results_brtp   = "/Users/danywaller/Projects/moon/equivalent_source_modeling/reiner_gamma_3body.txt"              # Fortran output
output_name  = "gravmag_sphere_map_2x2"

figsize = (10, 8)
save_png = True
DPI = 300

# Color scaling
symmetric_clim = True         # center colorbar at 0 using max(|min|,|max|)
override_clim = None          # override (vmin, vmax), e.g. (-50, 50)

# Plot style
use_pcolormesh = True         # True uses pcolormesh on inferred grid, False uses scatter

# ============================================================
# Helpers
# ============================================================

def _split_tokens(line: str):
    # .in should have no inline comments, but this is a safeguard
    for marker in ("[", "!", "#"):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.split()


def parse_input_cards(path: str):
    """
    Parse Cards 1–5 needed for plotting:
      Card 1: title
      Card 2: lat0 lon0 dlat dlon elvo nlat nlon
      Card 5: ifield (field selector)
    """
    with open(path, "r") as f:
        raw = [ln.strip() for ln in f if ln.strip() != ""]

    # Card 1: title (A80)
    title = raw[0]

    # Card 2
    t2 = _split_tokens(raw[1])
    if len(t2) < 7:
        raise ValueError(f"Card 2 must have 7 fields: lat lon dlat dlon elvo nlat nlon. Got: {raw[1]}")
    lat0, lon0, dlat, dlon, elvo = map(float, t2[:5])
    nlat, nlon = map(int, t2[5:7])

    # Card 5: ifield nrem rf1 ri1 rd1
    t5 = _split_tokens(raw[4])
    if len(t5) < 2:
        raise ValueError(f"Card 5 must start with: ifield nrem ... Got: {raw[4]}")
    ifield = int(float(t5[0]))  # tolerate "2.0" etc

    return title, lat0, lon0, dlat, dlon, elvo, nlat, nlon, ifield


def field_label_from_ifield(ifield: int):
    """
    Map ifield to label consistent with Fortran.

    ifield:
      0 -> gravity anomaly
      1 -> total magnetic field
      2 -> Br
      3 -> Bphi
      4 -> Btheta
    """
    if ifield == 0:
        return "gravity", "mGal"
    if ifield == 1:
        return "Btot", "nT"
    if ifield == 2:
        return "Br", "nT"
    if ifield == 3:
        return "Bphi", "nT"
    if ifield == 4:
        return "Btheta", "nT"
    return f"ifield={ifield}", "nT"


def infer_regular_grid(lon, lat, val):
    """
    Infer a regular (lat, lon) grid and reshape values accordingly.
    Returns: ulon, ulat, Z with Z shape (nlat, nlon), where rows correspond to ulat ascending.
    """
    ulon = np.unique(lon)
    ulat = np.unique(lat)

    nlon = ulon.size
    nlat = ulat.size

    if nlon * nlat != val.size:
        raise ValueError(
            "Cannot infer a complete regular grid: "
            "len(unique(lon))*len(unique(lat)) != Npoints. "
            "Switch to scatter (use_pcolormesh=False) or fix missing points."
        )

    lon_index = {x: i for i, x in enumerate(ulon)}
    lat_index = {y: i for i, y in enumerate(ulat)}

    Z = np.full((nlat, nlon), np.nan, dtype=float)
    for x, y, v in zip(lon, lat, val):
        Z[lat_index[y], lon_index[x]] = v

    if np.isnan(Z).any():
        raise ValueError("Grid has missing cells (NaNs).")

    return ulon, ulat, Z


def compute_clim(values):
    if override_clim is not None:
        return override_clim
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if symmetric_clim:
        m = max(abs(vmin), abs(vmax))
        return (-m, m)
    return (vmin, vmax)


def compute_clim_btot(values):
    """
    For |B| use [0, max] scaling
    """
    vmax = float(np.nanmax(values))
    return (0.0, vmax)


def load_results_multi(path: str):
    """
    Load either:
      (A) legacy 3-column: lon lat value
      (B) new 6-column:    lon lat Br Btheta Bphi Btot
    Skips comment lines starting with '#'.
    """
    xyz = np.loadtxt(path, comments="#")
    if xyz.ndim != 2 or xyz.shape[1] < 3:
        raise ValueError("results file must have at least 3 columns: lon lat value")

    if xyz.shape[1] >= 6:
        # new format :)
        lon = xyz[:, 0].astype(float)
        lat = xyz[:, 1].astype(float)
        br = xyz[:, 2].astype(float)
        btheta = xyz[:, 3].astype(float)
        bphi = xyz[:, 4].astype(float)
        btot = xyz[:, 5].astype(float)
        fields = {
            "Br": br,
            "Btheta": btheta,
            "Bphi": bphi,
            "|B|": btot,
        }
        return lon, lat, fields, "multi"
    else:
        # legacy format (single field)
        lon = xyz[:, 0].astype(float)
        lat = xyz[:, 1].astype(float)
        val = xyz[:, 2].astype(float)
        fields = {"value": val}
        return lon, lat, fields, "single"


def plot_panel(ax, lon, lat, values, lon0, lon1, lat0, lat1, title, *, cmap="RdBu_r", clim=None):
    if clim is None:
        clim = compute_clim(values)

    if use_pcolormesh:
        ulon, ulat, Z = infer_regular_grid(lon, lat, values)
        im = ax.pcolormesh(
            ulon, ulat, Z,
            shading="nearest",
            cmap=cmap,
            vmin=clim[0],
            vmax=clim[1],
        )
    else:
        im = ax.scatter(
            lon, lat,
            c=values,
            s=12,
            cmap=cmap,
            vmin=clim[0],
            vmax=clim[1],
            edgecolors="none",
        )

    ax.set_title(title)
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="box")
    return im


# ============================================================
# Parse input + load results
# ============================================================

title, lat0, lon0, dlat, dlon, elvo, nlat, nlon, ifield = parse_input_cards(input_in_file)
field_name, unit = field_label_from_ifield(ifield)

# Derived region (must match Fortran!)
lat1 = lat0 + (nlat - 1) * dlat
lon1 = lon0 + (nlon - 1) * dlon

print("Parsed from input file:")
print(f"  Title:  {title}")
print(f"  ifield: {ifield} -> {field_name} [{unit}]")
print(f"  Lat range: {lat0} → {lat1}")
print(f"  Lon range: {lon0} → {lon1}")
print(f"  Spacing:   {dlon}° × {dlat}°")
print(f"  Grid size: {nlat} × {nlon}")

lon, lat, fields, mode = load_results_multi(results_brtp)

print(f"Loaded results from {results_brtp} (mode={mode})")
print(f"  Npoints: {lon.size}")
print(f"  Lon range (results): [{lon.min():.3f}, {lon.max():.3f}]")
print(f"  Lat range (results): [{lat.min():.3f}, {lat.max():.3f}]")

# ============================================================
# Plot
# ============================================================

os.makedirs("./output/", exist_ok=True)

if mode == "multi":
    # 2x2: Br, Btheta, Bphi, |B|
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.ravel()

    panel_order = [
        ("Br", "$B_r$"),
        ("Btheta", "$B_θ$"),
        ("Bphi", "$B_φ$"),
        ("|B|", "|B|"),
    ]

    for ax, (key, label) in zip(axes, panel_order):

        if key == "|B|":
            # |B|: jet, 0 -> max
            clim = compute_clim_btot(fields[key])
            cmap = "jet"
        else:
            # components: diverging, symmetric around 0
            clim = None
            cmap = "RdBu_r"

        im = plot_panel(
            ax=ax,
            lon=lon,
            lat=lat,
            values=fields[key],
            lon0=lon0, lon1=lon1, lat0=lat0, lat1=lat1,
            title=label,
            cmap=cmap,
            clim=clim,
        )

        # colorbar beneath each subplot
        cbar = fig.colorbar(
            im,
            ax=ax,
            orientation="horizontal",
            pad=0.2,
            shrink=0.5,
        )
        cbar.set_label("nT")

    # Axis labels only on outer edges
    for ax in axes:
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")

    fig.suptitle(title, y=0.98)

    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_png = f"./output/{output_name}_Br_Btheta_Bphi_Btot.png"
    if save_png:
        fig.savefig(out_png, dpi=DPI)

    plt.show()
    print(f"Wrote: {out_png}")

else:
    # Legacy single-field behavior (1 panel)
    fig, ax = plt.subplots(figsize=(8, 4.8))

    im = plot_panel(
        ax=ax,
        lon=lon,
        lat=lat,
        values=fields["value"],
        lon0=lon0, lon1=lon1, lat0=lat0, lat1=lat1,
        title=f"{title} | {field_name}",
        unit=unit,
    )

    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")

    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label(unit)

    fig.tight_layout()

    out_png = f"./output/{output_name}_{field_name}.png"
    if save_png:
        fig.savefig(out_png, dpi=DPI)
    plt.show()
    print(f"Wrote: {out_png}")