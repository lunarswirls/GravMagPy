#!/usr/bin/env python3
"""plot summed cartesian magnetic components for one gravmag sphere case"""

import os
from pathlib import Path

import matplotlib

if os.environ.get("DISPLAY", "") == "" and os.environ.get("MPLBACKEND", "") == "":
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .utils.paths import artifact_path, observation_input_path


def read_noncomment_lines(path):
    """return stripped nonempty input lines excluding comments"""
    lines = []
    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if line and not line.startswith(("#", "!")):
                lines.append(line)
    return lines


def read_input(path):
    """read observation-grid metadata and all source outlines"""
    lines = read_noncomment_lines(path)
    bodies = []
    index = 0
    grid = None

    while index < len(lines):
        if index + 6 >= len(lines):
            raise ValueError(f"incomplete body block near line {index + 1} in {path}")

        title = lines[index]
        card2 = lines[index + 1].split()
        card3 = lines[index + 2].split()
        card5 = lines[index + 4].split()
        if len(card2) < 7 or len(card3) < 4 or len(card5) < 5:
            raise ValueError(f"invalid card block for {title!r} in {path}")

        body_grid = {
            "lat0": float(card2[0]),
            "lon0": float(card2[1]),
            "dlat": float(card2[2]),
            "dlon": float(card2[3]),
            "nlat": int(card2[5]),
            "nlon": int(card2[6]),
        }
        if grid is None:
            grid = body_grid
        elif body_grid != grid:
            raise ValueError("all bodies must use identical observation grids")

        nblim = int(card3[3])
        amplitude = float(card5[2])
        inclination = float(card5[3])
        geometry = lines[index + 6].split()

        if nblim == 1:
            if len(geometry) < 6:
                raise ValueError(f"invalid fixed-limit geometry for {title!r}")
            lat_max, lat_min, lon_max, lon_min = map(float, geometry[:4])
            lat = np.array([lat_min, lat_max, lat_max, lat_min, lat_min])
            lon = np.array([lon_min, lon_min, lon_max, lon_max, lon_min])
            index += 7
        else:
            if len(geometry) < 3:
                raise ValueError(f"invalid polygon header for {title!r}")
            npoints = int(geometry[0])
            vertices = lines[index + 7:index + 7 + npoints]
            if len(vertices) != npoints:
                raise ValueError(f"missing polygon vertices for {title!r}")
            coordinates = np.array([[float(value) for value in row.split()[:2]] for row in vertices])
            lat = coordinates[:, 0]
            lon = coordinates[:, 1]
            index += 7 + npoints

        bodies.append(
            {
                "title": title,
                "lon": lon,
                "lat": lat,
                "amplitude": amplitude,
                "inclination": inclination,
            }
        )

    if not bodies or grid is None:
        raise ValueError(f"no source bodies found in {path}")
    return {"grid": grid, "bodies": bodies}


def read_output_header(path):
    """return the concatenated output comment header"""
    headers = []
    with path.open("r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.strip()
            if line.startswith("#"):
                headers.append(line[1:].strip())
            elif line:
                break
    return " ".join(headers)


def read_and_sum_output(path):
    """load solver output and sum vector components over body ids"""
    header = read_output_header(path)
    if header and not all(component in header for component in ("Bx", "By", "Bz")):
        raise ValueError(f"expected cartesian Bx/By/Bz columns in {path}")

    table = np.loadtxt(path, comments="#")
    if table.ndim == 1:
        table = table[np.newaxis, :]
    if table.ndim != 2 or table.shape[1] < 7:
        raise ValueError(f"expected seven-column solver output in {path}")
    if not np.isfinite(table[:, :7]).all():
        raise ValueError(f"nonfinite values found in {path}")

    lon = np.round(table[:, 1], 6)
    lat = np.round(table[:, 2], 6)
    points = np.column_stack((lon, lat))
    unique_points, inverse = np.unique(points, axis=0, return_inverse=True)
    fields = np.zeros((unique_points.shape[0], 3), dtype=float)
    for column in range(3):
        np.add.at(fields[:, column], inverse, table[:, column + 3])

    total = np.sqrt(np.sum(fields * fields, axis=1))
    unit = "nT" if "_nT" in header else ""
    return {
        "lon": unique_points[:, 0],
        "lat": unique_points[:, 1],
        "bx": fields[:, 0],
        "by": fields[:, 1],
        "bz": fields[:, 2],
        "btot": total,
        "unit": unit,
    }


def make_grid(data, input_data):
    """reshape summed output into the regular grid declared by card 2"""
    lon_values = np.unique(data["lon"])
    lat_values = np.unique(data["lat"])
    lon_values.sort()
    lat_values.sort()

    expected_nlon = input_data["grid"]["nlon"]
    expected_nlat = input_data["grid"]["nlat"]
    if lon_values.size != expected_nlon or lat_values.size != expected_nlat:
        raise ValueError(
            "output grid dimensions do not match card 2: "
            f"found {lat_values.size} by {lon_values.size}, "
            f"expected {expected_nlat} by {expected_nlon}"
        )

    lon_index = {value: index for index, value in enumerate(lon_values)}
    lat_index = {value: index for index, value in enumerate(lat_values)}
    grids = {}
    for key in ("bx", "by", "bz", "btot"):
        grid = np.full((lat_values.size, lon_values.size), np.nan)
        for lon, lat, value in zip(data["lon"], data["lat"], data[key]):
            grid[lat_index[lat], lon_index[lon]] = value
        if np.isnan(grid).any():
            raise ValueError(f"output has missing grid cells for {key}")
        grids[key] = grid

    return lon_values, lat_values, grids


def add_source_outlines(axis, bodies):
    """draw source polygons over a field panel"""
    for body in bodies:
        axis.plot(body["lon"], body["lat"], color="white", linewidth=1.8, alpha=0.85)
        axis.plot(body["lon"], body["lat"], color="black", linewidth=0.7, alpha=0.9)


def plot_fields(input_path, output_path, image_path=None):
    """render the summed Bx, By, Bz, and Btot maps"""
    input_path, output_path = map(Path, (input_path, output_path))
    image_path = (artifact_path(input_path, kind="figs", name=f"{output_path.stem}_bxyz.png")
                  if image_path is None else Path(image_path).expanduser().resolve())
    input_data = read_input(input_path)
    data = read_and_sum_output(output_path)
    lon, lat, grids = make_grid(data, input_data)

    figure, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    labels = {"bx": "Bx", "by": "By", "bz": "Bz", "btot": "Btot"}

    for axis, key in zip(axes.ravel(), ("bx", "by", "bz", "btot")):
        grid = grids[key]
        if key == "btot":
            lower = 0.0
            upper = float(np.nanmax(grid))
            upper = upper if upper > 0.0 else 1.0
            colormap = "viridis"
        else:
            limit = float(np.nanmax(np.abs(grid)))
            limit = limit if limit > 0.0 else 1.0
            lower = -limit
            upper = limit
            colormap = "RdBu_r"

        image = axis.pcolormesh(
            lon,
            lat,
            grid,
            shading="nearest",
            cmap=colormap,
            vmin=lower,
            vmax=upper,
        )
        add_source_outlines(axis, input_data["bodies"])
        axis.set_title(labels[key])
        axis.set_xlabel("longitude [deg]")
        axis.set_ylabel("latitude [deg]")
        axis.set_aspect("equal", adjustable="box")
        colorbar = figure.colorbar(image, ax=axis, orientation="horizontal", pad=0.12, shrink=0.85)
        colorbar.set_label(data["unit"] or "field units")

    case_name = input_path.stem.replace("_", " ")
    body_count = len(input_data["bodies"])
    body_label = "body" if body_count == 1 else "bodies"
    figure.suptitle(f"{case_name}\nsummed response from {body_count} source {body_label}")
    image_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(image_path, dpi=250)
    plt.close(figure)
    return image_path


def plot_fit(observations, result, image_path=None):
    """plot measured and modeled orbital components against sample index"""
    if image_path is None:
        image_path = artifact_path(observation_input_path(observations), kind="figs", suffix="_fit.png")
    observed = np.asarray(observations["field_nt"])
    predicted = np.asarray(result["predicted_nt"])
    if predicted.shape != observed.shape:
        raise ValueError("predictions and observations must have identical shapes")
    figure, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True, constrained_layout=True)
    for index, (axis, label) in enumerate(zip(axes, ("Bx", "By", "Bz", "Btot"))):
        actual = observed[:, index] if index < 3 else np.linalg.norm(observed, axis=1)
        model = predicted[:, index] if index < 3 else np.linalg.norm(predicted, axis=1)
        axis.plot(actual, color="0.25", linewidth=0.7, label="observed")
        axis.plot(model, color="tab:orange", linewidth=0.8, label="modeled")
        axis.set_ylabel(f"{label} [nT]")
        axis.grid(alpha=0.15)
    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel("observation index (input file and row order)")
    figure.suptitle(f"orbital source fit | rmse={result['rmse_nt']:.3g} nT | success={result['success']}")
    image_path = Path(image_path).expanduser().resolve()
    image_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(image_path, dpi=180)
    plt.close(figure)
    return image_path
