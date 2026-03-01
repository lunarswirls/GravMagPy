#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# edit me! :)
# ============================================================
input_in_file = "gravmag_sphere_1body_test.in"
results_brtp  = "gravmag_sphere_1body_test.txt"
output_name   = "gravmag_sphere_map_2x2"

figsize = (10, 8)
save_png = True
DPI = 300

# Color scaling
symmetric_clim = True
override_clim = None

# Plot style
use_pcolormesh = True

# Overlay bodies (outlines from .in file)
overlay_bodies = True
label_bodies = False
outline_lw = 1.5

# ============================================================
# choose plotting mode for multi-body outputs
# ============================================================
plot_each_body_separately = False
plot_total_sum = True
# ============================================================

# ============================================================
# Helpers
# ============================================================

def _split_tokens(line: str):
    for marker in ("[", "!", "#"):
        if marker in line:
            line = line.split(marker, 1)[0]
    return line.split()


def wrap180(lon_deg: float) -> float:
    """Wrap to [-180, 180)."""
    x = (lon_deg + 180.0) % 360.0 - 180.0
    # keep +180 exclusive
    if x >= 180.0:
        x -= 360.0
    return x


def unwrap_lons_to_ref(lons: np.ndarray, lon_ref: float) -> np.ndarray:
    """
    Shift longitudes by +/-360 so they lie on a continuous branch near lon_ref.
    Inputs/outputs in degrees; lon_ref assumed already in [-180,180).
    """
    out = np.array([wrap180(x) for x in lons], dtype=float)
    for i in range(out.size):
        while out[i] - lon_ref >= 180.0:
            out[i] -= 360.0
        while out[i] - lon_ref < -180.0:
            out[i] += 360.0
    return out


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
    vmax = float(np.nanmax(values))
    return (0.0, vmax)


def _detect_results_layout(xyz: np.ndarray):
    """
    Supports:
      old single:  lon lat value                -> (N,3)
      old multi :  lon lat Br Btheta Bphi Btot  -> (N,6)
      new single:  body lon lat value           -> (N,4)
      new multi :  body lon lat Bx By Bz Btot   -> (N,7+)
    """
    if xyz.ndim != 2:
        raise ValueError("results file must be 2D")

    ncol = xyz.shape[1]
    if ncol == 3:
        return False, "single", dict(lon=0, lat=1, val0=2)
    if ncol == 6:
        return False, "multi", dict(lon=0, lat=1, br=2, bth=3, bph=4, btot=5)
    if ncol == 4:
        return True, "single", dict(body=0, lon=1, lat=2, val0=3)
    if ncol >= 7:
        return True, "multi", dict(body=0, lon=1, lat=2, bx=3, by=4, bz=5, btot=6)

    raise ValueError(f"Unrecognized results column count: {ncol}")


def load_results(path: str):
    xyz = np.loadtxt(path, comments="#")
    if xyz.ndim != 2 or xyz.shape[1] < 3:
        raise ValueError("results file must have at least 3 numeric columns")

    has_body, mode, idx = _detect_results_layout(xyz)

    body = None
    if has_body:
        body = xyz[:, idx["body"]].astype(int)

    lon = xyz[:, idx["lon"]].astype(float)
    lat = xyz[:, idx["lat"]].astype(float)

    if mode == "multi":
        # Old 6-col mode (Br/Btheta/Bphi/Btot) still supported
        if "bx" in idx:
            fields = {
                "Bx": xyz[:, idx["bx"]].astype(float),
                "By": xyz[:, idx["by"]].astype(float),
                "Bz": xyz[:, idx["bz"]].astype(float),
                "|B|": xyz[:, idx["btot"]].astype(float),
            }
        else:
            fields = {
                "Br": xyz[:, idx["br"]].astype(float),
                "Btheta": xyz[:, idx["bth"]].astype(float),
                "Bphi": xyz[:, idx["bph"]].astype(float),
                "|B|": xyz[:, idx["btot"]].astype(float),
            }
    else:
        fields = {"value": xyz[:, idx["val0"]].astype(float)}

    return body, lon, lat, fields, mode


def bin_to_regular_grid(lon, lat, val, *, decimals=6):
    """
    Robust gridding:
      - round lon/lat to a fixed decimal
      - build full grid
      - average duplicates
      - error if missing cells after binning
    Returns sorted (ulon, ulat, Z) where Z is (nlat,nlon).
    """
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)
    val = np.asarray(val, float)

    lon_r = np.round(lon, decimals=decimals)
    lat_r = np.round(lat, decimals=decimals)

    ulon = np.unique(lon_r)
    ulat = np.unique(lat_r)
    ulon.sort()
    ulat.sort()

    nlon = ulon.size
    nlat = ulat.size

    lon_index = {x: i for i, x in enumerate(ulon)}
    lat_index = {y: i for i, y in enumerate(ulat)}

    Zsum = np.zeros((nlat, nlon), dtype=float)
    Zcnt = np.zeros((nlat, nlon), dtype=int)

    for x, y, v in zip(lon_r, lat_r, val):
        ii = lat_index[y]
        jj = lon_index[x]
        Zsum[ii, jj] += v
        Zcnt[ii, jj] += 1

    Z = np.full((nlat, nlon), np.nan, dtype=float)
    m = Zcnt > 0
    Z[m] = Zsum[m] / Zcnt[m]

    if np.isnan(Z).any():
        missing = np.count_nonzero(np.isnan(Z))
        raise ValueError(f"Grid has missing cells after binning (missing={missing}). "
                         f"Try lowering decimals (currently {decimals}).")

    return ulon, ulat, Z


def plot_panel(ax, lon, lat, values, title, *, cmap="RdBu_r", clim=None, grid_decimals=6):
    if clim is None:
        clim = compute_clim(values)

    if use_pcolormesh:
        ulon, ulat, Z = bin_to_regular_grid(lon, lat, values, decimals=grid_decimals)
        im = ax.pcolormesh(
            ulon, ulat, Z,
            shading="nearest",
            cmap=cmap,
            vmin=clim[0],
            vmax=clim[1],
        )
        # Set limits from grid axes (most reliable)
        ax.set_xlim(float(ulon.min()), float(ulon.max()))
        ax.set_ylim(float(ulat.min()), float(ulat.max()))
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
        ax.set_xlim(float(np.nanmin(lon)), float(np.nanmax(lon)))
        ax.set_ylim(float(np.nanmin(lat)), float(np.nanmax(lat)))

    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="box")
    return im


# ============================================================
# Option-B .in parser (Card 7 polygon header always)
# ============================================================

def parse_all_bodies(path: str):
    """
    Parses repeated body blocks until EOF.

    Option-B Card 7:
      npts depth_top_km depth_bot_km
      lat lon   (repeated npts)

    Returns:
      bodies: list of dicts with polygon outline in lon/lat degrees (wrapped [-180,180))
    """
    with open(path, "r") as f:
        raw = [ln.strip() for ln in f if ln.strip()]

    bodies = []
    idx = 0

    while idx < len(raw):
        title = raw[idx]; idx += 1

        # Cards 2-6
        if idx + 5 > len(raw):
            raise ValueError(f"Unexpected EOF while reading cards 2-6 for body {len(bodies)+1}")

        t2 = _split_tokens(raw[idx]); idx += 1
        lat0, lon0, dlat, dlon, elvo = map(float, t2[:5])
        nlat, nlon = map(int, t2[5:7])

        t3 = _split_tokens(raw[idx]); idx += 1
        nr, ntheta, nphi, nblim = map(int, t3[:4])

        t4 = _split_tokens(raw[idx]); idx += 1
        htheta, hphi, phi1, phi2, rho = map(float, t4[:5])

        t5 = _split_tokens(raw[idx]); idx += 1
        ifield = int(float(t5[0]))
        nrem   = int(float(t5[1]))
        mamp   = float(t5[2]) if len(t5) > 2 else 0.0
        inc    = float(t5[3]) if len(t5) > 3 else 0.0
        dec    = float(t5[4]) if len(t5) > 4 else 0.0

        t6 = _split_tokens(raw[idx]); idx += 1
        iprint = int(float(t6[0]))
        nfile  = int(float(t6[1]))
        nopt   = int(float(t6[2]))
        first  = float(t6[3]) if len(t6) > 3 else 0.0
        conint = float(t6[4]) if len(t6) > 4 else 0.0

        # Card 7 (Option B): polygon header + vertices
        if idx >= len(raw):
            raise ValueError(f"Unexpected EOF while reading Card 7 for body {len(bodies)+1}")

        t7 = _split_tokens(raw[idx]); idx += 1
        if len(t7) < 3:
            raise ValueError(f"Card 7 must be: npts depth_top_km depth_bot_km. Got: {raw[idx-1]}")

        npts = int(float(t7[0]))
        depth_top_km = float(t7[1])
        depth_bot_km = float(t7[2])

        if npts < 3:
            raise ValueError(f"Polygon needs >=3 vertices, got npts={npts}")

        if idx + npts > len(raw):
            raise ValueError(f"Unexpected EOF while reading {npts} polygon vertices for body {len(bodies)+1}")

        plat = np.zeros(npts, dtype=float)
        plon = np.zeros(npts, dtype=float)

        for k in range(npts):
            tp = _split_tokens(raw[idx]); idx += 1
            if len(tp) < 2:
                raise ValueError(f"Polygon vertex must be: lat lon. Got: {raw[idx-1]}")
            plat[k] = float(tp[0])
            plon[k] = wrap180(float(tp[1]))

        # unwrap lon branch like Fortran (mean lon_ref)
        lon_ref = wrap180(float(np.mean(plon)))
        plon_u = unwrap_lons_to_ref(plon, lon_ref)

        # close polygon if needed
        if not (np.isclose(plon_u[0], plon_u[-1]) and np.isclose(plat[0], plat[-1])):
            plon_u = np.r_[plon_u, plon_u[0]]
            plat   = np.r_[plat,   plat[0]]

        bodies.append(
            dict(
                title=title,
                lat0=lat0, lon0=lon0, dlat=dlat, dlon=dlon, elvo=elvo, nlat=nlat, nlon=nlon,
                nr=nr, ntheta=ntheta, nphi=nphi, nblim=nblim,
                htheta=htheta, hphi=hphi, phi1=phi1, phi2=phi2, rho=rho,
                ifield=ifield, nrem=nrem, M_amp=mamp, inc=inc, dec=dec,
                iprint=iprint, nfile=nfile, nopt=nopt, first=first, conint=conint,
                depth_top_km=depth_top_km, depth_bot_km=depth_bot_km,
                poly_lon=plon_u, poly_lat=plat,
            )
        )

    return bodies


def assert_common_grid(bodies):
    b0 = bodies[0]
    keys = ("lat0","lon0","dlat","dlon","elvo","nlat","nlon")
    for bi, b in enumerate(bodies[1:], start=2):
        for k in keys:
            if float(b[k]) != float(b0[k]):
                print(f"WARNING: Body {bi} differs in Card 2 field '{k}': {b[k]} vs {b0[k]}")
    return b0


def overlay_body_outlines(ax, bodies):
    for bi, b in enumerate(bodies, start=1):
        ax.plot(b["poly_lon"], b["poly_lat"], color="magenta", linewidth=outline_lw, label=f"Body {bi}")
        if label_bodies and b["poly_lon"].size > 0:
            cx = float(np.mean(b["poly_lon"]))
            cy = float(np.mean(b["poly_lat"]))
            ax.text(
                cx, cy, f"{bi}",
                ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
            )


def plot_multi_2x2(lon, lat, fields, title, bodies, out_png, *, grid_decimals=6):
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.ravel()

    panel_order = [
        ("Bx", "$B_x$"),
        ("By", "$B_y$"),
        ("Bz", "$B_z$"),
        ("|B|", "|B|"),
    ]

    for ax, (key, label) in zip(axes, panel_order):
        if key == "|B|":
            clim = compute_clim_btot(fields[key])
            cmap = "jet"
        else:
            clim = None
            cmap = "RdBu_r"

        im = plot_panel(ax, lon, lat, fields[key], label, cmap=cmap, clim=clim, grid_decimals=grid_decimals)

        if overlay_bodies:
            overlay_body_outlines(ax, bodies)

        cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.2, shrink=0.5)
        cbar.set_label("nT")

        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")

    fig.suptitle(title, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    if save_png:
        fig.savefig(out_png, dpi=DPI)
    plt.show()
    print(f"Wrote: {out_png}")


def plot_single_panel(lon, lat, values, title, bodies, unit, out_png, *, grid_decimals=6):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    im = plot_panel(ax, lon, lat, values, title, cmap="RdBu_r", clim=None, grid_decimals=grid_decimals)

    if overlay_bodies:
        overlay_body_outlines(ax, bodies)

    ax.set_xlabel("Longitude [deg]")
    ax.set_ylabel("Latitude [deg]")
    cbar = fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label(unit)

    fig.tight_layout()
    if save_png:
        fig.savefig(out_png, dpi=DPI)
    plt.show()
    print(f"Wrote: {out_png}")


def aggregate_sum_over_bodies(body_id, lon, lat, fields, *, decimals=6):
    lon = np.asarray(lon, float)
    lat = np.asarray(lat, float)

    lon_r = np.round(lon, decimals=decimals)
    lat_r = np.round(lat, decimals=decimals)

    pts = np.column_stack((lon_r, lat_r))
    uniq_pts, inv = np.unique(pts, axis=0, return_inverse=True)

    lon_u = uniq_pts[:, 0].astype(float)
    lat_u = uniq_pts[:, 1].astype(float)

    comp_keys = [k for k in ("Bx", "By", "Bz") if k in fields]
    out = {}
    for k in comp_keys:
        v = np.asarray(fields[k], float)
        acc = np.zeros(uniq_pts.shape[0], dtype=float)
        np.add.at(acc, inv, v)
        out[k] = acc

    if all(k in out for k in ("Bx", "By", "Bz")):
        out["|B|"] = np.sqrt(out["Bx"]**2 + out["By"]**2 + out["Bz"]**2)

    return lon_u, lat_u, out


# ============================================================
# Main
# ============================================================

bodies = parse_all_bodies(input_in_file)
if len(bodies) == 0:
    raise RuntimeError("No bodies parsed from input file.")

b0 = assert_common_grid(bodies)

body_id, lon, lat, fields, mode = load_results(results_brtp)

print("Parsed bodies from input file:")
print(f"  Nbodies: {len(bodies)}")
for i, b in enumerate(bodies, start=1):
    print(f"  Body {i}: title='{b['title']}'  nblim={b['nblim']}  depth={b['depth_top_km']}..{b['depth_bot_km']} km")

print(f"Loaded results from {results_brtp} (mode={mode})")
print(f"  Npoints: {lon.size}")
print(f"  body_id column present: {body_id is not None}")

os.makedirs("./output/", exist_ok=True)

if (not plot_each_body_separately) and (not plot_total_sum):
    plot_total_sum = True

# Prefer plot extents from DATA (robust to wrap/rounding)
lon_min, lon_max = float(np.nanmin(lon)), float(np.nanmax(lon))
lat_min, lat_max = float(np.nanmin(lat)), float(np.nanmax(lat))
print(f"Data extents: lon {lon_min}..{lon_max}  lat {lat_min}..{lat_max}")

if mode == "multi":
    # per-body plots
    if plot_each_body_separately:
        if body_id is None:
            raise RuntimeError("plot_each_body_separately=True requires results file with body_id column.")
        for bid in np.unique(body_id):
            m = (body_id == bid)
            title = f"Body {bid}: {bodies[bid-1]['title']}" if (1 <= bid <= len(bodies)) else f"Body {bid}"
            out_png = f"./output/{output_name}_body{bid:02d}_Bx_By_Bz_Btot.png"
            plot_multi_2x2(
                lon[m], lat[m],
                {k: v[m] for k, v in fields.items()},
                title, bodies, out_png,
                grid_decimals=6,
            )

    # total/summed plot
    if plot_total_sum:
        if body_id is not None and len(np.unique(body_id)) > 1:
            lon_s, lat_s, fields_s = aggregate_sum_over_bodies(body_id, lon, lat, fields, decimals=6)
            title = f"Total (sum of {len(np.unique(body_id))} bodies): {b0['title']}"
            out_png = f"./output/{output_name}_TOTALSUM_Bx_By_Bz_Btot.png"
            plot_multi_2x2(lon_s, lat_s, fields_s, title, bodies, out_png, grid_decimals=6)
        else:
            title = f"Total: {b0['title']}"
            out_png = f"./output/{output_name}_TOTAL_Bx_By_Bz_Btot.png"
            plot_multi_2x2(lon, lat, fields, title, bodies, out_png, grid_decimals=6)

else:
    # single-field mode
    unit = "nT"
    if plot_each_body_separately:
        if body_id is None:
            raise RuntimeError("plot_each_body_separately=True requires results file with body_id column.")
        for bid in np.unique(body_id):
            m = (body_id == bid)
            title = f"Body {bid}: {bodies[bid-1]['title']} | value" if (1 <= bid <= len(bodies)) else f"Body {bid} | value"
            out_png = f"./output/{output_name}_body{bid:02d}_value.png"
            plot_single_panel(lon[m], lat[m], fields["value"][m], title, bodies, unit, out_png, grid_decimals=6)

    if plot_total_sum:
        title = f"Total: {b0['title']} | value"
        out_png = f"./output/{output_name}_TOTAL_value.png"
        plot_single_panel(lon, lat, fields["value"], title, bodies, unit, out_png, grid_decimals=6)