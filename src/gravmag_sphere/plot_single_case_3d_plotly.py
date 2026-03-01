#!/usr/bin/env python3
"""
Create a 3D Plotly visualization for a single source case.

Inputs:
- solver input file (.in) for one body (first body block is used)
- BRTP output file with columns: body_id lon lat Br Btheta Bphi Btot (or without body_id)

Output:
- self-contained HTML with:
  - source geometry wireframe (top, bottom, side edges)
  - sparse Br/Btheta/Bphi projected vectors in 3D
  - vector colors mapped by Btot
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def read_noncomment_lines(path: Path) -> List[str]:
    """Return non-empty, non-comment lines from an input-card file."""
    lines: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("!"):
                continue
            lines.append(s)
    return lines


def parse_single_source_input(path: Path) -> Dict[str, object]:
    """Parse the first body block from one .in file for geometry metadata."""
    lines = read_noncomment_lines(path)
    if len(lines) < 7:
        raise RuntimeError(f"Could not parse first body cards from {path}")

    title = lines[0]
    c2 = lines[1].split()
    c3 = lines[2].split()
    c7 = lines[6].split()

    lat0_deg = float(c2[0])
    lon0_deg = float(c2[1])
    dlat_deg = float(c2[2])
    dlon_deg = float(c2[3])
    elvo_km = float(c2[4])
    nlat = int(float(c2[5]))
    nlon = int(float(c2[6]))

    nblim = int(float(c3[3]))
    depth_top_km: float
    depth_bot_km: float
    poly_lat: List[float]
    poly_lon: List[float]

    if nblim == 1:
        # fixed limits: lat_max lat_min lon_max lon_min depth_top depth_bot
        if len(c7) < 6:
            raise RuntimeError("Invalid Card 7 fixed-limits line.")
        lat_max = float(c7[0])
        lat_min = float(c7[1])
        lon_max = float(c7[2])
        lon_min = float(c7[3])
        depth_top_km = float(c7[4])
        depth_bot_km = float(c7[5])
        poly_lat = [lat_min, lat_max, lat_max, lat_min]
        poly_lon = [lon_min, lon_min, lon_max, lon_max]
    else:
        # polygon: npts depth_top depth_bot + npts vertices
        if len(c7) < 3:
            raise RuntimeError("Invalid Card 7 polygon header.")
        npts = int(float(c7[0]))
        depth_top_km = float(c7[1])
        depth_bot_km = float(c7[2])
        if len(lines) < 7 + npts:
            raise RuntimeError(f"Polygon expects {npts} vertices, but input is too short.")
        poly_lat = []
        poly_lon = []
        for i in range(npts):
            vals = lines[7 + i].split()
            if len(vals) < 2:
                raise RuntimeError(f"Invalid polygon vertex at line index {7+i}.")
            p1 = float(vals[0])
            p2 = float(vals[1])
            # Accept both (lat lon) and (lon lat)
            if abs(p1) <= 90.0 and abs(p2) > 90.0:
                lat, lon = p1, p2
            elif abs(p1) > 90.0 and abs(p2) <= 90.0:
                lon, lat = p1, p2
            else:
                lat, lon = p1, p2
            poly_lat.append(lat)
            poly_lon.append(lon)

    return {
        "title": title,
        "lat0_deg": lat0_deg,
        "lon0_deg": lon0_deg,
        "dlat_deg": dlat_deg,
        "dlon_deg": dlon_deg,
        "elvo_km": elvo_km,
        "nlat": nlat,
        "nlon": nlon,
        "depth_top_km": depth_top_km,
        "depth_bot_km": depth_bot_km,
        "poly_lat": np.array(poly_lat, dtype=float),
        "poly_lon": np.array(poly_lon, dtype=float),
        "nblim": nblim,
    }


def parse_brtp(path: Path) -> Dict[str, np.ndarray]:
    """Load BRTP rows with flexible support for with/without body_id."""
    unit = ""
    lon: List[float] = []
    lat: List[float] = []
    br: List[float] = []
    bt: List[float] = []
    bp: List[float] = []
    btot: List[float] = []
    body: List[int] = []

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "_nT" in s:
                    unit = "nT"
                elif "_mGal" in s or "_MGAL" in s:
                    unit = "mGal"
                continue

            vals = s.split()
            if len(vals) < 6:
                continue
            try:
                if len(vals) >= 7:
                    body_id = int(float(vals[0]))
                    xlon = float(vals[1])
                    xlat = float(vals[2])
                    xbr = float(vals[3])
                    xbt = float(vals[4])
                    xbp = float(vals[5])
                    xmag = float(vals[6])
                else:
                    body_id = 0
                    xlon = float(vals[0])
                    xlat = float(vals[1])
                    xbr = float(vals[2])
                    xbt = float(vals[3])
                    xbp = float(vals[4])
                    xmag = float(vals[5])
            except ValueError:
                continue

            body.append(body_id)
            lon.append(xlon)
            lat.append(xlat)
            br.append(xbr)
            bt.append(xbt)
            bp.append(xbp)
            btot.append(xmag)

    if not lon:
        raise RuntimeError(f"No data rows found in {path}")

    return {
        "unit": unit,
        "body_id": np.array(body, dtype=int),
        "lon": np.array(lon, dtype=float),
        "lat": np.array(lat, dtype=float),
        "Br": np.array(br, dtype=float),
        "Btheta": np.array(bt, dtype=float),
        "Bphi": np.array(bp, dtype=float),
        "Btot": np.array(btot, dtype=float),
    }


def lonlat_to_xyz_km(r_km: float, lat_deg: np.ndarray, lon_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert lon/lat coordinates on a sphere to cartesian xyz in km."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = r_km * np.cos(lat) * np.cos(lon)
    y = r_km * np.cos(lat) * np.sin(lon)
    z = r_km * np.sin(lat)
    return x, y, z


def basis_vectors(lat_deg: np.ndarray, lon_deg: np.ndarray) -> Dict[str, np.ndarray]:
    """Build local spherical basis vectors (er, etheta, ephi) in cartesian form."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    c_lat = np.cos(lat)
    s_lat = np.sin(lat)
    c_lon = np.cos(lon)
    s_lon = np.sin(lon)

    # r outward
    er = np.column_stack((c_lat * c_lon, c_lat * s_lon, s_lat))
    # theta southward (increasing colatitude)
    et = np.column_stack((s_lat * c_lon, s_lat * s_lon, -c_lat))
    # phi eastward
    ep = np.column_stack((-s_lon, c_lon, np.zeros_like(s_lon)))
    return {"Br": er, "Btheta": et, "Bphi": ep}


def sparse_indices(lon: np.ndarray, lat: np.ndarray, step: int, max_vectors: int) -> np.ndarray:
    """Pick a sparse but geographically distributed subset of grid points."""
    lon_r = np.round(lon, 6)
    lat_r = np.round(lat, 6)
    ulon = np.unique(lon_r)
    ulat = np.unique(lat_r)
    sel_lon = set(ulon[:: max(1, step)])
    sel_lat = set(ulat[:: max(1, step)])
    mask = np.array([(x in sel_lon) and (y in sel_lat) for x, y in zip(lon_r, lat_r)], dtype=bool)
    idx = np.where(mask)[0]
    if idx.size == 0:
        idx = np.arange(0, lon.size, max(1, step))
    if idx.size > max_vectors:
        pick = np.linspace(0, idx.size - 1, max_vectors).astype(int)
        idx = idx[pick]
    return idx


def ensure_closed(lat: np.ndarray, lon: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Close polygon arrays by repeating the first vertex at the end if needed."""
    if lat.size == 0:
        return lat, lon
    if (lat[0] != lat[-1]) or (lon[0] != lon[-1]):
        lat = np.append(lat, lat[0])
        lon = np.append(lon, lon[0])
    return lat, lon


def add_sphere_wireframe(fig, radius_km: float, go):
    """Add a lightweight reference sphere wireframe to the 3-D figure."""
    lat_lines = np.linspace(-60.0, 60.0, 7)
    lon_lines = np.linspace(-180.0, 180.0, 181)
    for lat0 in lat_lines:
        lat = np.full_like(lon_lines, lat0)
        x, y, z = lonlat_to_xyz_km(radius_km, lat, lon_lines)
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color="rgba(120,120,120,0.25)", width=1),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    lon_lines2 = np.linspace(-150.0, 150.0, 7)
    lat_lines2 = np.linspace(-90.0, 90.0, 181)
    for lon0 in lon_lines2:
        lon = np.full_like(lat_lines2, lon0)
        x, y, z = lonlat_to_xyz_km(radius_km, lat_lines2, lon)
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color="rgba(120,120,120,0.25)", width=1),
                showlegend=False,
                hoverinfo="skip",
            )
        )


def build_geometry_traces(fig, meta: Dict[str, object], rsphere_km: float, go):
    """Plot top, bottom, and side-edge traces for the source body geometry."""
    poly_lat = np.array(meta["poly_lat"], dtype=float)
    poly_lon = np.array(meta["poly_lon"], dtype=float)
    poly_lat, poly_lon = ensure_closed(poly_lat, poly_lon)

    r_top = rsphere_km - float(meta["depth_top_km"])
    r_bot = rsphere_km - float(meta["depth_bot_km"])

    xt, yt, zt = lonlat_to_xyz_km(r_top, poly_lat, poly_lon)
    xb, yb, zb = lonlat_to_xyz_km(r_bot, poly_lat, poly_lon)

    fig.add_trace(
        go.Scatter3d(
            x=xt,
            y=yt,
            z=zt,
            mode="lines",
            line=dict(color="#ff8c00", width=6),
            name="Source top boundary",
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=xb,
            y=yb,
            z=zb,
            mode="lines",
            line=dict(color="#cc6f00", width=5),
            name="Source bottom boundary",
        )
    )

    # side edges
    for i in range(max(0, len(xt) - 1)):
        fig.add_trace(
            go.Scatter3d(
                x=[xt[i], xb[i]],
                y=[yt[i], yb[i]],
                z=[zt[i], zb[i]],
                mode="lines",
                line=dict(color="rgba(220,120,0,0.8)", width=3),
                showlegend=False,
                hoverinfo="skip",
            )
        )


def make_human_title(input_title: str, output_path: Path) -> str:
    """Compose a readable figure title from input title + output stem."""
    stem = output_path.stem
    if stem.endswith("_brtp"):
        stem = stem[:-5]
    return f"{input_title} | {stem}"


def main() -> int:
    """CLI entrypoint for one-case 3-D Plotly visualization."""
    ap = argparse.ArgumentParser(
        description="Create a 3D Plotly HTML for one source case with sparse Br/Btheta/Bphi vector projections."
    )
    ap.add_argument("input_file", help="Single-case input .in file")
    ap.add_argument("output_file", help="Single-case BRTP output file")
    ap.add_argument("html_file", help="Output HTML path")
    ap.add_argument("--rsphere-km", type=float, default=1737.4)
    ap.add_argument("--sparse-step", type=int, default=10, help="Grid decimation step for vectors")
    ap.add_argument("--max-vectors", type=int, default=140, help="Maximum sparse vectors")
    ap.add_argument(
        "--vector-scale",
        type=float,
        default=0.03,
        help="Max component projection length as fraction of sphere radius",
    )
    args = ap.parse_args()

    try:
        import plotly.graph_objects as go
        from plotly.colors import sample_colorscale
    except Exception:
        print(
            "ERROR: Plotly is required for this script. Install it with:\n"
            "  python -m pip install plotly",
            file=sys.stderr,
        )
        return 2

    in_path = Path(args.input_file)
    out_path = Path(args.output_file)
    html_path = Path(args.html_file)
    if html_path.parent == Path("."):
        html_path = Path("figs") / html_path

    if not in_path.exists():
        raise FileNotFoundError(in_path)
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    html_path.parent.mkdir(parents=True, exist_ok=True)

    meta = parse_single_source_input(in_path)
    data = parse_brtp(out_path)

    unique_bodies = np.unique(data["body_id"])
    if unique_bodies.size > 1:
        print(
            "WARNING: output contains multiple body IDs; plotting all rows as provided.",
            file=sys.stderr,
        )

    ro_km = float(args.rsphere_km) + float(meta["elvo_km"])
    x0, y0, z0 = lonlat_to_xyz_km(ro_km, data["lat"], data["lon"])
    start = np.column_stack((x0, y0, z0))
    basis = basis_vectors(data["lat"], data["lon"])

    idx = sparse_indices(data["lon"], data["lat"], args.sparse_step, args.max_vectors)
    br = data["Br"][idx]
    bt = data["Btheta"][idx]
    bp = data["Bphi"][idx]
    btot = data["Btot"][idx]

    ref = max(
        float(np.max(np.abs(br))) if br.size else 0.0,
        float(np.max(np.abs(bt))) if bt.size else 0.0,
        float(np.max(np.abs(bp))) if bp.size else 0.0,
        1.0e-30,
    )
    comp_scale = float(args.vector_scale) * ro_km / ref

    vmin = float(np.nanmin(btot))
    vmax = float(np.nanmax(btot))
    if math.isclose(vmin, vmax):
        vmin -= 1.0
        vmax += 1.0
    norm = np.clip((btot - vmin) / (vmax - vmin), 0.0, 1.0)
    colors = [sample_colorscale("Viridis", float(t))[0] for t in norm]

    fig = go.Figure()
    add_sphere_wireframe(fig, ro_km, go)
    build_geometry_traces(fig, meta, float(args.rsphere_km), go)

    components = [("Br", br, 6), ("Btheta", bt, 5), ("Bphi", bp, 4)]
    legend_shown = {k: False for k, _, _ in components}

    for i_local, i_global in enumerate(idx):
        p0 = start[i_global]
        btot_i = btot[i_local]
        color_i = colors[i_local]
        for name, vals, width in components:
            vec = basis[name][i_global] * (vals[i_local] * comp_scale)
            p1 = p0 + vec
            hover = (
                f"{name}: {vals[i_local]:.4g}<br>"
                f"Btot: {btot_i:.4g} {data['unit'] or 'units'}<br>"
                f"lon: {data['lon'][i_global]:.3f} lat: {data['lat'][i_global]:.3f}"
            )
            fig.add_trace(
                go.Scatter3d(
                    x=[p0[0], p1[0]],
                    y=[p0[1], p1[1]],
                    z=[p0[2], p1[2]],
                    mode="lines",
                    line=dict(color=color_i, width=width),
                    name=f"{name} projection",
                    showlegend=not legend_shown[name],
                    hovertemplate=hover + "<extra></extra>",
                )
            )
            legend_shown[name] = True

    # Invisible marker field to carry colorbar for Btot
    fig.add_trace(
        go.Scatter3d(
            x=start[idx, 0],
            y=start[idx, 1],
            z=start[idx, 2],
            mode="markers",
            marker=dict(
                size=2,
                color=btot,
                colorscale="Viridis",
                cmin=vmin,
                cmax=vmax,
                showscale=True,
                colorbar=dict(title=f"Btot [{data['unit'] or 'units'}]"),
                opacity=0.01,
            ),
            hoverinfo="skip",
            showlegend=False,
            name="Btot scale",
        )
    )

    fig.update_layout(
        title=make_human_title(str(meta["title"]), out_path),
        scene=dict(
            xaxis_title="X [km]",
            yaxis_title="Y [km]",
            zaxis_title="Z [km]",
            aspectmode="data",
        ),
        legend=dict(x=0.01, y=0.99),
        margin=dict(l=0, r=0, b=0, t=45),
    )

    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(html_path), include_plotlyjs=True)
    print(f"Wrote 3D HTML: {html_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
