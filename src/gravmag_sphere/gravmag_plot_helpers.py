#!/usr/bin/env python3
"""
Utilities for loading, documenting, and plotting grav/mag Br/Btheta/Bphi outputs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import re

import os
import matplotlib
if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_brtp(path: Path) -> Dict[str, np.ndarray]:
    """Load BRTP table into named numpy arrays."""
    arr = np.loadtxt(path, comments="#")
    if arr.ndim != 2 or arr.shape[1] < 7:
        raise ValueError(f"Unexpected BRTP format in {path}")
    return {
        "body_id": arr[:, 0].astype(int),
        "lon": arr[:, 1],
        "lat": arr[:, 2],
        "Br": arr[:, 3],
        "Btheta": arr[:, 4],
        "Bphi": arr[:, 5],
        "Btot": arr[:, 6],
    }


def pretty_case_name(stem: str) -> str:
    """Convert compact case stem to a human-readable phrase."""
    s = stem
    if s.startswith("gravmag_sphere_"):
        s = s[len("gravmag_sphere_") :]

    s = re.sub(r"(\d+)body", r"\1 body", s)
    s = s.replace("fixedlim", "fixed limits")
    s = s.replace("complexlarge", "large complex polygon")
    s = s.replace("incmix", "mixed inc")
    s = s.replace("decmix", "mixed dec")
    s = re.sub(r"inc(-?\d+)", r"inc \1°", s)
    s = re.sub(r"dec(-?\d+)", r"dec \1°", s)

    token_map = {
        "mag": "magnetic",
        "grav": "gravity",
        "orient": "orientation",
        "base": "baseline",
        "weak": "weak",
        "polygon": "polygon",
        "fixed": "fixed",
        "incna": "inc n/a",
        "decna": "dec n/a",
    }

    toks = []
    for tok in s.replace("_", " ").split():
        toks.append(token_map.get(tok, tok))

    out = " ".join(toks)
    out = re.sub(r"\s+", " ", out).strip()
    if out:
        out = out[0].upper() + out[1:]
    return out


def parse_unit_from_header(path: Path) -> str:
    """Parse the first header block and infer display units."""
    unit = ""
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "_nT" in s:
                    unit = "nT"
                elif "_mGal" in s or "_MGAL" in s:
                    unit = "mGal"
                continue
            break
    return unit


def aggregate_total(data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Sum component vectors over body_id at each lon/lat sample."""
    lon_r = np.round(data["lon"], 6)
    lat_r = np.round(data["lat"], 6)
    pts = np.column_stack((lon_r, lat_r))
    uniq, inv = np.unique(pts, axis=0, return_inverse=True)

    br = np.zeros(uniq.shape[0], dtype=float)
    bt = np.zeros(uniq.shape[0], dtype=float)
    bp = np.zeros(uniq.shape[0], dtype=float)

    np.add.at(br, inv, data["Br"])
    np.add.at(bt, inv, data["Btheta"])
    np.add.at(bp, inv, data["Bphi"])
    btot = np.sqrt(br * br + bt * bt + bp * bp)

    return {"lon": uniq[:, 0], "lat": uniq[:, 1], "Br": br, "Btheta": bt, "Bphi": bp, "Btot": btot}


def to_grid(lon: np.ndarray, lat: np.ndarray, values: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Regrid 1-D lon/lat/value arrays into a dense 2-D map."""
    lon_r = np.round(lon, 6)
    lat_r = np.round(lat, 6)
    ulon = np.unique(lon_r)
    ulat = np.unique(lat_r)
    ulon.sort()
    ulat.sort()

    i_lon = {x: i for i, x in enumerate(ulon)}
    i_lat = {y: i for i, y in enumerate(ulat)}

    z = np.full((ulat.size, ulon.size), np.nan, dtype=float)
    for x, y, v in zip(lon_r, lat_r, values):
        z[i_lat[y], i_lon[x]] = v

    if np.isnan(z).any():
        raise RuntimeError(f"Grid contains {int(np.isnan(z).sum())} missing cells.")
    return ulon, ulat, z


def plot_brtp_2x2(
    data: Dict[str, np.ndarray], title: str, unit: str, *, figsize: Tuple[float, float] = (11, 8)
) -> plt.Figure:
    """Render Br/Btheta/Bphi/Btot in a 2x2 panel."""
    keys = ["Br", "Btheta", "Bphi", "Btot"]
    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)
    axes = axes.ravel()

    for ax, key in zip(axes, keys):
        ulon, ulat, z = to_grid(data["lon"], data["lat"], data[key])
        if key == "Btot":
            vmin, vmax = 0.0, float(np.nanmax(z))
            cmap = "viridis"
        else:
            m = float(np.nanmax(np.abs(z)))
            vmin, vmax = -m, m
            cmap = "RdBu_r"

        im = ax.pcolormesh(ulon, ulat, z, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(key)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_aspect("equal", adjustable="box")
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.12, shrink=0.85)
        cb.set_label(unit if unit else "units")

    fig.suptitle(title)
    return fig


def _noncomment_lines(path: Path) -> List[str]:
    """Read non-empty, non-comment lines from an input-card file."""
    lines = []
    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            s = raw.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith("!"):
                continue
            lines.append(s)
    return lines


def describe_example(input_path: Path) -> Dict[str, str]:
    """Summarize the first body block from an example input file."""
    lines = _noncomment_lines(input_path)
    if len(lines) < 7:
        raise RuntimeError(f"Could not parse cards from {input_path}")

    title = lines[0]
    c2 = lines[1].split()
    c3 = lines[2].split()
    c5 = lines[4].split()
    ifield = int(float(c5[0]))
    mode = "magnetic" if ifield == 2 else "gravity"

    grid = (
        f"obs grid: lat0={c2[0]} lon0={c2[1]} dlat={c2[2]} dlon={c2[3]} "
        f"nlat={c2[5]} nlon={c2[6]}"
    )
    mesh = f"card3 mesh seed: nr={c3[0]} ntheta={c3[1]} nphi={c3[2]} nblim={c3[3]}"

    prop = (
        f"magnetization={c5[2]} A/m, inc={c5[3]} deg, dec={c5[4]} deg"
        if ifield == 2
        else f"density contrast={c5[2]} kg/m^3"
    )

    nblim = int(float(c3[3]))
    if nblim == 1:
        g = lines[6].split()
        geom = (
            "fixed limits: "
            f"lat[{g[1]},{g[0]}], lon[{g[3]},{g[2]}], depth_top={g[4]} km, depth_bot={g[5]} km"
        )
    else:
        head = lines[6].split()
        npts = int(float(head[0]))
        geom = f"polygon: npts={npts}, depth_top={head[1]} km, depth_bot={head[2]} km"

    return {"title": title, "mode": mode, "property": prop, "geometry": geom, "grid": grid, "mesh": mesh}


def plot_example(
    example_stem: str,
    *,
    root_dir: Path | str = ".",
    examples_dir: str = "examples",
    output_dir: str = "output",
) -> plt.Figure:
    """Convenience wrapper: load one example BRTP file and return its 2x2 figure."""
    root = Path(root_dir)
    input_path = root / examples_dir / f"{example_stem}.in"
    brtp_path = root / output_dir / f"{example_stem}_brtp.txt"
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    if not brtp_path.exists():
        raise FileNotFoundError(brtp_path)

    meta = describe_example(input_path)
    data = load_brtp(brtp_path)
    total = aggregate_total(data)
    unit = parse_unit_from_header(brtp_path)
    title = f"{pretty_case_name(example_stem)}: Br / Btheta / Bphi / Btot ({meta['mode']})"
    return plot_brtp_2x2(total, title, unit)
