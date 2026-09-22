#!/usr/bin/env python3
from pathlib import Path
import os

import numpy as np
from scipy.interpolate import RegularGridInterpolator


HERE = Path(__file__).resolve().parent
CACHE_DIR = HERE / ".plot_cache"
MPL_DIR = CACHE_DIR / "matplotlib"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(MPL_DIR)

CASE_STEM = "reiner_gamma_table4_surface_direct"
INPUT_IN = HERE / f"{CASE_STEM}.in"
INPUT_BRTP = HERE / "output" / f"{CASE_STEM}_brtp.txt"
OUTPUT_PNG = HERE / "fig" / f"{CASE_STEM}_topography_3d.png"
OUTPUT_HTML = HERE / "fig" / f"{CASE_STEM}_topography_3d.html"
TOPO_CACHE_NPZ = HERE / "output" / f"{CASE_STEM}_topography_patch.npz"

MOON_MEAN_RADIUS_KM = 1737.4
TOPO_LMAX = 512
SURFACE_STEP = 2
VERTICAL_EXAGGERATION = 22.0
HILLSHADE_EXAGGERATION = 38.0
VIEW_ELEVATION_DEG = 32.0
VIEW_AZIMUTH_DEG = -128.0
FIGURE_SIZE = (18.0, 15.5)
PLOTLY_WIDTH_PX = 1700
PLOTLY_HEIGHT_PX = 1350
SOURCE_MARKER_SIZE = 18.0
SOURCE_OFFSET_KM = 0.02
LIGHT_AZIMUTH_DEG = 315.0
LIGHT_ALTITUDE_DEG = 42.0


def read_noncomment_lines(path: Path) -> list[str]:
    lines: list[str] = []
    with open(path, "r", encoding="utf-8") as stream:
        for raw in stream:
            text = raw.strip()
            if not text or text.startswith("#") or text.startswith("!"):
                continue
            lines.append(text)
    return lines


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


def parse_first_grid_meta(path: Path) -> dict[str, float]:
    lines = read_noncomment_lines(path)
    if len(lines) < 2:
        raise RuntimeError(f"Could not parse grid metadata from {path}")
    card2 = lines[1].split()
    return {
        "lat0_deg": float(card2[0]),
        "lon0_deg": float(card2[1]),
        "dlat_deg": float(card2[2]),
        "dlon_deg": float(card2[3]),
        "elvo_km": float(card2[4]),
        "nlat": int(float(card2[5])),
        "nlon": int(float(card2[6])),
    }


def parse_fixed_box_centers(path: Path) -> np.ndarray:
    lines = read_noncomment_lines(path)
    centers: list[tuple[float, float]] = []
    idx = 0
    while idx < len(lines):
        if idx + 6 >= len(lines):
            raise RuntimeError(f"Incomplete body block near line {idx + 1} in {path}")
        card3 = lines[idx + 2].split()
        nblim = int(float(card3[3]))
        if nblim != 1:
            raise RuntimeError("This script expects fixed-limit bodies in the Reiner Gamma input.")
        card7 = lines[idx + 6].split()
        if len(card7) < 6:
            raise RuntimeError(f"Invalid Card 7 fixed-limit line near line {idx + 7}")
        lat_max = float(card7[0])
        lat_min = float(card7[1])
        lon_max = float(card7[2])
        lon_min = float(card7[3])
        centers.append((0.5 * (lat_min + lat_max), 0.5 * (lon_min + lon_max)))
        idx += 7
    return np.asarray(centers, dtype=float)


def wrap_lon_360(lon_deg: np.ndarray) -> np.ndarray:
    return np.mod(lon_deg, 360.0)


def build_topography_patch(lat_grid_deg: np.ndarray, lon_grid_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        import pyshtools as pysh
    except Exception as exc:  # pragma: no cover - environment-specific import
        raise RuntimeError(
            "pyshtools is required to build the topography cache. "
            "Run this script once with /Users/danywaller/code/GravMagPy/.venv_compare/bin/python."
        ) from exc

    coeffs = pysh.datasets.Moon.LDEM_shape_pa()
    topo_grid = coeffs.expand(lmax=TOPO_LMAX)

    nlat = topo_grid.nlat
    nlon = topo_grid.nlon
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon
    lat_axis_deg = (90.0 - np.arange(nlat) * dlat)[::-1]
    lon_axis_deg = np.arange(nlon) * dlon
    radius_km = np.flipud(topo_grid.data) / 1000.0

    lon_axis_ext = np.concatenate([lon_axis_deg, [360.0]])
    radius_ext = np.concatenate([radius_km, radius_km[:, :1]], axis=1)

    topo_interp = RegularGridInterpolator(
        (lat_axis_deg, lon_axis_ext),
        radius_ext,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    lon_grid_360 = wrap_lon_360(lon_grid_deg)
    radius_patch_km = topo_interp(
        np.column_stack([lat_grid_deg.ravel(), lon_grid_360.ravel()])
    ).reshape(lat_grid_deg.shape)
    elev_patch_km = radius_patch_km - MOON_MEAN_RADIUS_KM

    np.savez_compressed(
        TOPO_CACHE_NPZ,
        lat=lat_grid_deg,
        lon=lon_grid_deg,
        radius_km=radius_patch_km,
        elev_km=elev_patch_km,
    )
    return radius_patch_km, elev_patch_km


def load_or_build_topography_patch(
    lat_grid_deg: np.ndarray,
    lon_grid_deg: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, str]:
    if TOPO_CACHE_NPZ.is_file():
        cache = np.load(TOPO_CACHE_NPZ)
        lat_cache = cache["lat"]
        lon_cache = cache["lon"]
        if lat_cache.shape == lat_grid_deg.shape and lon_cache.shape == lon_grid_deg.shape:
            if np.allclose(lat_cache, lat_grid_deg) and np.allclose(lon_cache, lon_grid_deg):
                return cache["radius_km"], cache["elev_km"], "cache"
    radius_km, elev_km = build_topography_patch(lat_grid_deg, lon_grid_deg)
    return radius_km, elev_km, "pyshtools"


def source_plot_positions(
    centers: np.ndarray,
    lat0_deg: float,
    lon0_deg: float,
    topo_elev_interp: RegularGridInterpolator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    surf_elev_km = topo_elev_interp(centers)
    _ = lat0_deg
    _ = lon0_deg
    lon_deg = centers[:, 1]
    lat_deg = centers[:, 0]
    z_km = surf_elev_km * VERTICAL_EXAGGERATION + SOURCE_OFFSET_KM
    return lon_deg, lat_deg, z_km


def field_limits_and_cmap(key: str, data: np.ndarray) -> tuple[float, float, str]:
    if key == "Btot":
        return 0.0, float(np.nanmax(data)), "Viridis"
    vmax = float(np.nanmax(np.abs(data)))
    return -vmax, vmax, "RdBu_r"


def write_static_png(
    fields: dict[str, np.ndarray],
    lon_ds: np.ndarray,
    lat_ds: np.ndarray,
    z_ds: np.ndarray,
    elev_ds: np.ndarray,
    src_lon_deg: np.ndarray,
    src_lat_deg: np.ndarray,
    src_z_km: np.ndarray,
    altitude_km: float,
    topo_source_label: str,
) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        from matplotlib import cm, colors
        from matplotlib.colors import LightSource
        import matplotlib.pyplot as plt
    except Exception:
        return False

    def field_norm_and_cmap(key: str, data: np.ndarray) -> tuple[colors.Normalize, object]:
        if key == "Btot":
            vmax = float(np.nanmax(data))
            return colors.Normalize(vmin=0.0, vmax=vmax), cm.viridis
        vmax = float(np.nanmax(np.abs(data)))
        return colors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax), cm.RdBu_r

    def style_axis(ax, title: str) -> None:
        ax.set_title(title, pad=12.0)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")
        ax.set_zlabel(f"Elevation x{VERTICAL_EXAGGERATION:.0f} [km]")
        ax.view_init(elev=VIEW_ELEVATION_DEG, azim=VIEW_AZIMUTH_DEG)
        ax.set_xlim(float(np.nanmin(lon_ds)), float(np.nanmax(lon_ds)))
        ax.set_ylim(float(np.nanmin(lat_ds)), float(np.nanmax(lat_ds)))
        ax.set_zlim(float(np.nanmin(z_ds)) - 0.4, float(np.nanmax(z_ds)) + 0.4)
        ax.set_box_aspect(
            (
                float(np.nanmax(lon_ds) - np.nanmin(lon_ds)),
                float(np.nanmax(lat_ds) - np.nanmin(lat_ds)),
                max(float(np.nanmax(z_ds) - np.nanmin(z_ds)), 1.0),
            )
        )

    ls = LightSource(azdeg=LIGHT_AZIMUTH_DEG, altdeg=LIGHT_ALTITUDE_DEG)
    keys = ["Br", "Btheta", "Bphi", "Btot"]
    fig = plt.figure(figsize=FIGURE_SIZE)

    for panel_idx, key in enumerate(keys, start=1):
        ax = fig.add_subplot(2, 2, panel_idx, projection="3d")
        data_ds = fields[key][::SURFACE_STEP, ::SURFACE_STEP]
        norm, cmap = field_norm_and_cmap(key, data_ds)
        rgb = cmap(norm(data_ds))[..., :3]
        shaded_rgb = ls.shade_rgb(rgb, elev_ds, vert_exag=HILLSHADE_EXAGGERATION, blend_mode="soft")

        ax.plot_surface(
            lon_ds,
            lat_ds,
            z_ds,
            rstride=1,
            cstride=1,
            facecolors=shaded_rgb,
            linewidth=0.0,
            antialiased=False,
            shade=False,
        )
        ax.scatter(
            src_lon_deg,
            src_lat_deg,
            src_z_km,
            s=SOURCE_MARKER_SIZE,
            c="black",
            edgecolors="white",
            linewidths=0.5,
            alpha=0.95,
            depthshade=False,
        )

        title = (
            f"{key}\n"
            f"min {np.nanmin(fields[key]):.1f} nT, max {np.nanmax(fields[key]):.1f} nT"
        )
        style_axis(ax, title)

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, shrink=0.62, pad=0.06, label=f"{key} [nT]")

    fig.suptitle(
        "Reiner Gamma Table 4 forward field on lunar topography\n"
        f"{CASE_STEM}, altitude {altitude_km:.1f} km, "
        f"{src_lon_deg.size} discrete sources, topography {topo_source_label} (lmax {TOPO_LMAX})",
        fontsize=16,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return True


def write_plotly_html(
    fields: dict[str, np.ndarray],
    lon_ds: np.ndarray,
    lat_ds: np.ndarray,
    z_ds: np.ndarray,
    src_lon_deg: np.ndarray,
    src_lat_deg: np.ndarray,
    src_z_km: np.ndarray,
    altitude_km: float,
    topo_source_label: str,
) -> bool:
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except Exception:
        return False

    panel_specs = [("Br", 1, 1), ("Btheta", 1, 2), ("Bphi", 2, 1), ("Btot", 2, 2)]
    fig = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"type": "scene"}, {"type": "scene"}], [{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=[
            f"{key}: {np.nanmin(fields[key]):.1f} to {np.nanmax(fields[key]):.1f} nT"
            for key, _, _ in panel_specs
        ],
        horizontal_spacing=0.02,
        vertical_spacing=0.03,
    )

    for key, row, col in panel_specs:
        data_ds = fields[key][::SURFACE_STEP, ::SURFACE_STEP]
        cmin, cmax, colorscale = field_limits_and_cmap(key, fields[key])
        if colorscale == "RdBu_r":
            colorscale = "RdBu"
            reversescale = True
        else:
            reversescale = False

        surface = go.Surface(
            x=lon_ds,
            y=lat_ds,
            z=z_ds,
            surfacecolor=data_ds,
            colorscale=colorscale,
            reversescale=reversescale,
            cmin=cmin,
            cmax=cmax,
            showscale=True,
            colorbar=dict(
                title=f"{key} [nT]",
                len=0.34,
                thickness=16,
            ),
            hovertemplate=(
                "Longitude %{x:.3f} deg<br>"
                "Latitude %{y:.3f} deg<br>"
                "Elev*exag %{z:.2f} km<br>"
                f"{key} %{{surfacecolor:.2f}} nT<extra></extra>"
            ),
            lighting=dict(ambient=0.75, diffuse=0.7, roughness=0.95, specular=0.06),
            lightposition=dict(x=350, y=250, z=900),
        )
        fig.add_trace(surface, row=row, col=col)
        fig.add_trace(
            go.Scatter3d(
                x=src_lon_deg,
                y=src_lat_deg,
                z=src_z_km,
                mode="markers",
                marker=dict(size=3, color="black", line=dict(color="white", width=1)),
                name="Source centers",
                showlegend=(row == 1 and col == 1),
                hovertemplate="Source center<extra></extra>",
            ),
            row=row,
            col=col,
        )

    x_span = float(np.nanmax(lon_ds) - np.nanmin(lon_ds))
    y_span = float(np.nanmax(lat_ds) - np.nanmin(lat_ds))
    z_span = max(float(np.nanmax(z_ds) - np.nanmin(z_ds)), 1.0)
    aspect_ratio = dict(x=x_span, y=y_span, z=z_span)
    camera = dict(eye=dict(x=1.45, y=-1.55, z=0.82))

    for scene_name in ["scene", "scene2", "scene3", "scene4"]:
        fig.layout[scene_name].update(
            xaxis_title="Longitude [deg]",
            yaxis_title="Latitude [deg]",
            zaxis_title=f"Elevation x{VERTICAL_EXAGGERATION:.0f} [km]",
            aspectmode="manual",
            aspectratio=aspect_ratio,
            camera=camera,
        )

    fig.update_layout(
        title=(
            "Reiner Gamma Table 4 forward field on lunar topography<br>"
            f"{CASE_STEM}, altitude {altitude_km:.1f} km, "
            f"{src_lon_deg.size} discrete sources, topography {topo_source_label} (lmax {TOPO_LMAX})"
        ),
        width=PLOTLY_WIDTH_PX,
        height=PLOTLY_HEIGHT_PX,
        margin=dict(l=0, r=0, b=0, t=80),
        legend=dict(x=0.01, y=0.99),
    )
    OUTPUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(OUTPUT_HTML), include_plotlyjs=True)
    return True


if not INPUT_IN.is_file():
    raise FileNotFoundError(f"Missing input file: {INPUT_IN}")

if not INPUT_BRTP.is_file():
    raise FileNotFoundError(
        f"Missing BRTP output: {INPUT_BRTP}. Run run_reiner_gamma_table4_surface_direct.sh first."
    )

grid_meta = parse_first_grid_meta(INPUT_IN)
fields = load_brtp_sum(INPUT_BRTP)
source_centers = parse_fixed_box_centers(INPUT_IN)

lon_grid_deg, lat_grid_deg = np.meshgrid(fields["lon"], fields["lat"])
radius_km, elev_km, topo_source_label = load_or_build_topography_patch(lat_grid_deg, lon_grid_deg)
z_plot_km = elev_km * VERTICAL_EXAGGERATION

topo_elev_interp = RegularGridInterpolator(
    (fields["lat"], fields["lon"]),
    elev_km,
    method="linear",
    bounds_error=False,
    fill_value=None,
)
lat0_deg = float(0.5 * (np.nanmin(fields["lat"]) + np.nanmax(fields["lat"])))
lon0_deg = float(0.5 * (np.nanmin(fields["lon"]) + np.nanmax(fields["lon"])))
src_lon_deg, src_lat_deg, src_z_km = source_plot_positions(
    source_centers,
    lat0_deg,
    lon0_deg,
    topo_elev_interp,
)

lon_ds = lon_grid_deg[::SURFACE_STEP, ::SURFACE_STEP]
lat_ds = lat_grid_deg[::SURFACE_STEP, ::SURFACE_STEP]
z_ds = z_plot_km[::SURFACE_STEP, ::SURFACE_STEP]
elev_ds = elev_km[::SURFACE_STEP, ::SURFACE_STEP]

png_written = write_static_png(
    fields,
    lon_ds,
    lat_ds,
    z_ds,
    elev_ds,
    src_lon_deg,
    src_lat_deg,
    src_z_km,
    grid_meta["elvo_km"],
    topo_source_label,
)
html_written = write_plotly_html(
    fields,
    lon_ds,
    lat_ds,
    z_ds,
    src_lon_deg,
    src_lat_deg,
    src_z_km,
    grid_meta["elvo_km"],
    topo_source_label,
)

print(f"Input IN: {INPUT_IN}")
print(f"Input BRTP: {INPUT_BRTP}")
print(f"Topography cache: {TOPO_CACHE_NPZ} ({topo_source_label})")
if png_written:
    print(f"Wrote 3D PNG: {OUTPUT_PNG}")
else:
    print("Skipped PNG: matplotlib unavailable.")
if html_written:
    print(f"Wrote 3D HTML: {OUTPUT_HTML}")
else:
    print("Skipped HTML: plotly unavailable.")
for key in ["Br", "Btheta", "Bphi", "Btot"]:
    print(
        f"{key}: min {np.nanmin(fields[key]):.9f} nT, "
        f"max {np.nanmax(fields[key]):.9f} nT"
    )
