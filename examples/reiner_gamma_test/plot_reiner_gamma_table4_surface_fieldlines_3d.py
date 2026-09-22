#!/usr/bin/env python3
from pathlib import Path
import math
import os

import numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator


HERE = Path(__file__).resolve().parent
CACHE_DIR = HERE / ".plot_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(CACHE_DIR / "matplotlib")

CASE_STEM = "reiner_gamma_table4_surface_direct"
INPUT_IN = HERE / f"{CASE_STEM}.in"
INPUT_BRTP = HERE / "output" / f"{CASE_STEM}_brtp.txt"
TOPO_CACHE_NPZ = HERE / "output" / f"{CASE_STEM}_topography_patch.npz"
WAC_CACHE_NPZ = HERE / "output" / f"{CASE_STEM}_wac_patch.npz"
WAC_MOSAIC_NAME = "Lunar_LRO_LROC-WAC_Mosaic_global_100m_June2013.tif"
WAC_MOSAIC_PATHS = (
    HERE / WAC_MOSAIC_NAME,
    Path("/Users/danywaller/Downloads") / WAC_MOSAIC_NAME,
)
OUTPUT_HTML = HERE / "figs" / f"{CASE_STEM}_fieldlines_3d.html"

MOON_MEAN_RADIUS_KM = 1737.4
MU0_OVER_4PI = 1.0e-7
PLOT_LON_MIN_DEG = -60.5
PLOT_LON_MAX_DEG = -56.5
PLOT_LAT_MIN_DEG = 5.5
PLOT_LAT_MAX_DEG = 9.0
PLOT_STEP_DEG = 0.02
SURFACE_DECIMATE = 1
ISO_BTOT_NT = 50.0
ALTITUDE_LEVELS_KM = tuple(np.arange(0.0, 30.0 + 1.0, 1.0))
ISO_MIN_COMPONENT_PIXELS = 12
ISO_SMOOTH_SIGMA_PX = 1.4
ISO_MASK_KEEP_THRESHOLD = 0.08
ISO_EDGE_TAPER_START = 0.20
SOURCE_SOFTENING_KM = 0.1
PLOTLY_WIDTH_PX = 800
PLOTLY_HEIGHT_PX = 600
SURFACE_OPACITY = 1.0
ISO_SURFACE_OPACITY = 0.6
SOURCE_MARKER_SIZE = 3


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

    def reshape(values: np.ndarray) -> np.ndarray:
        grid = np.full((len(lat_u), len(lon_u)), np.nan)
        lon_to_ix = {v: i for i, v in enumerate(lon_u)}
        lat_to_iy = {v: i for i, v in enumerate(lat_u)}
        for (lo, la), val in zip(uniq, values):
            grid[lat_to_iy[la], lon_to_ix[lo]] = val
        return grid

    return {
        "lon": lon_u,
        "lat": lat_u,
        "Br": reshape(br),
        "Btheta": reshape(bt),
        "Bphi": reshape(bp),
        "Btot": reshape(btot),
    }


def wrap_lon_360(lon_deg: np.ndarray) -> np.ndarray:
    return np.mod(lon_deg, 360.0)


def wrap_lon_180(lon_deg: np.ndarray) -> np.ndarray:
    return ((lon_deg + 180.0) % 360.0) - 180.0


def match_lon_domain(lon_deg: np.ndarray, ref_lon_deg: np.ndarray) -> np.ndarray:
    if float(np.nanmin(ref_lon_deg)) < 0.0:
        return wrap_lon_180(lon_deg)
    return wrap_lon_360(lon_deg)


def lonlatr_to_xyz_m(lat_deg: np.ndarray, lon_deg: np.ndarray, radius_m: np.ndarray) -> np.ndarray:
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    cos_lat = np.cos(lat_rad)
    return np.column_stack(
        (
            radius_m * cos_lat * np.cos(lon_rad),
            radius_m * cos_lat * np.sin(lon_rad),
            radius_m * np.sin(lat_rad),
        )
    )


def xyz_to_lonlatr_m(xyz_m: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radius_m = np.linalg.norm(xyz_m, axis=1)
    lon_deg = np.degrees(np.arctan2(xyz_m[:, 1], xyz_m[:, 0]))
    lat_deg = np.degrees(np.arcsin(np.clip(xyz_m[:, 2] / radius_m, -1.0, 1.0)))
    return lon_deg, lat_deg, radius_m


def local_basis(lat0_deg: float, lon0_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat = math.radians(lat0_deg)
    lon = math.radians(lon0_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0], dtype=float)
    north = np.array(
        [-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)],
        dtype=float,
    )
    up = np.array(
        [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)],
        dtype=float,
    )
    return east, north, up


def global_to_local_km(
    xyz_m: np.ndarray,
    origin_xyz_m: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    delta_m = xyz_m - origin_xyz_m.reshape(1, 3)
    x_km = delta_m @ east / 1000.0
    y_km = delta_m @ north / 1000.0
    z_km = delta_m @ up / 1000.0
    return x_km, y_km, z_km


def angle_distance_deg(lon_a_deg: float, lon_b_deg: float) -> float:
    return ((lon_a_deg - lon_b_deg + 180.0) % 360.0) - 180.0


def load_topography_patch(
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

    try:
        import pyshtools as pysh
    except Exception as exc:  # pragma: no cover - environment-specific import
        raise RuntimeError(
            "Expanded topographic basemap requires a matching cached patch or pyshtools. "
            "Build the cache once with /Users/danywaller/code/GravMagPy/.venv_compare/bin/python."
        ) from exc

    coeffs = pysh.datasets.Moon.LDEM_shape_pa()
    lmax = min(512, int(coeffs.lmax))
    topo_grid = coeffs.expand(lmax=lmax)
    nlat = topo_grid.nlat
    nlon = topo_grid.nlon
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon
    lat_axis_deg = (90.0 - np.arange(nlat) * dlat)[::-1]
    lon_axis_deg = np.arange(nlon) * dlon
    radius_grid_km = np.flipud(topo_grid.data) / 1000.0
    lon_axis_ext = np.concatenate([lon_axis_deg, [360.0]])
    radius_ext = np.concatenate([radius_grid_km, radius_grid_km[:, :1]], axis=1)
    interp = RegularGridInterpolator(
        (lat_axis_deg, lon_axis_ext),
        radius_ext,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )
    lon_360 = wrap_lon_360(lon_grid_deg)
    radius_km = interp(np.column_stack([lat_grid_deg.ravel(), lon_360.ravel()])).reshape(lat_grid_deg.shape)
    elev_km = radius_km - MOON_MEAN_RADIUS_KM
    np.savez_compressed(
        TOPO_CACHE_NPZ,
        lat=lat_grid_deg,
        lon=lon_grid_deg,
        radius_km=radius_km,
        elev_km=elev_km,
    )
    return radius_km, elev_km, "pyshtools"


def normalize_grid(grid: np.ndarray) -> np.ndarray:
    finite = np.isfinite(grid)
    if not np.any(finite):
        raise ValueError("Grid has no finite values to normalize.")
    lo = float(np.nanmin(grid))
    hi = float(np.nanmax(grid))
    if hi <= lo:
        return np.where(finite, 0.5, np.nan)
    return (grid - lo) / (hi - lo)


def load_wac_patch(
    lat_grid_deg: np.ndarray,
    lon_grid_deg: np.ndarray,
) -> tuple[np.ndarray, str]:
    if WAC_CACHE_NPZ.is_file():
        cache = np.load(WAC_CACHE_NPZ)
        lat_cache = cache["lat"]
        lon_cache = cache["lon"]
        if lat_cache.shape == lat_grid_deg.shape and lon_cache.shape == lon_grid_deg.shape:
            if np.allclose(lat_cache, lat_grid_deg) and np.allclose(lon_cache, lon_grid_deg):
                return cache["wac_norm"], "cache"

    wac_path = None
    for candidate in WAC_MOSAIC_PATHS:
        if candidate.is_file():
            wac_path = candidate
            break
    if wac_path is None:
        raise FileNotFoundError(
            "Missing WAC mosaic. Checked: " + ", ".join(str(path) for path in WAC_MOSAIC_PATHS)
        )

    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.warp import Resampling, reproject

    try:
        src = rasterio.open(wac_path)
    except Exception as exc:
        raise RuntimeError(
            f"Unable to read WAC mosaic at {wac_path}. "
            "Grant Codex access to Downloads or place a copy next to this script in the workspace."
        ) from exc

    with src:
        if src.crs is None:
            raise ValueError("WAC mosaic GeoTIFF has no CRS.")

        dst_crs = src.crs.geodetic_crs or src.crs
        dst_transform = from_bounds(
            float(np.nanmin(lon_grid_deg)),
            float(np.nanmin(lat_grid_deg)),
            float(np.nanmax(lon_grid_deg)),
            float(np.nanmax(lat_grid_deg)),
            lon_grid_deg.shape[1],
            lat_grid_deg.shape[0],
        )

        if src.count == 1:
            patch = np.full(lat_grid_deg.shape, np.nan, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=patch,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=dst_crs,
                dst_nodata=np.nan,
                resampling=Resampling.bilinear,
            )
            wac = np.flipud(patch)
        else:
            bands = []
            for band_index in range(1, min(src.count, 3) + 1):
                patch = np.full(lat_grid_deg.shape, np.nan, dtype=np.float32)
                reproject(
                    source=rasterio.band(src, band_index),
                    destination=patch,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    dst_nodata=np.nan,
                    resampling=Resampling.bilinear,
                )
                bands.append(np.flipud(patch))
            while len(bands) < 3:
                bands.append(bands[-1].copy())
            wac = 0.2126 * bands[0] + 0.7152 * bands[1] + 0.0722 * bands[2]

    wac_norm = normalize_grid(wac)
    np.savez_compressed(WAC_CACHE_NPZ, lat=lat_grid_deg, lon=lon_grid_deg, wac_norm=wac_norm)
    return wac_norm, f"rasterio:{wac_path}"


def parse_body_models(path: Path) -> list[dict[str, float]]:
    lines = read_noncomment_lines(path)
    bodies: list[dict[str, float]] = []
    idx = 0
    while idx < len(lines):
        if idx + 6 >= len(lines):
            raise RuntimeError(f"Incomplete body block near line {idx + 1} in {path}")
        title = lines[idx]
        _ = title
        card3 = lines[idx + 2].split()
        card5 = lines[idx + 4].split()
        card7 = lines[idx + 6].split()
        nblim = int(float(card3[3]))
        if nblim != 1:
            raise RuntimeError("This script expects fixed-limit Reiner Gamma sources.")
        if len(card5) < 5 or len(card7) < 6:
            raise RuntimeError(f"Invalid body definition near line {idx + 1} in {path}")

        mag_amp_apm = float(card5[2])
        inc_deg = float(card5[3])
        dec_deg = float(card5[4])

        lat_max = float(card7[0])
        lat_min = float(card7[1])
        lon_max = float(card7[2])
        lon_min = float(card7[3])
        depth_top_km = float(card7[4])
        depth_bot_km = float(card7[5])

        lon_c_deg = 0.5 * (lon_min + lon_max)
        lat_c_deg = 0.5 * (lat_min + lat_max)

        lon_width_rad = math.radians(lon_max - lon_min)
        lat_band = math.sin(math.radians(lat_max)) - math.sin(math.radians(lat_min))
        r_top_m = (MOON_MEAN_RADIUS_KM - depth_top_km) * 1000.0
        r_bot_m = (MOON_MEAN_RADIUS_KM - depth_bot_km) * 1000.0
        volume_m3 = abs(lon_width_rad * lat_band * (r_top_m**3 - r_bot_m**3) / 3.0)
        r_centroid_m = 0.75 * abs(r_top_m**4 - r_bot_m**4) / abs(r_top_m**3 - r_bot_m**3)

        cos_i = math.cos(math.radians(inc_deg))
        sin_i = math.sin(math.radians(inc_deg))
        cos_d = math.cos(math.radians(dec_deg))
        sin_d = math.sin(math.radians(dec_deg))
        moment_unit = np.array([cos_i * cos_d, cos_i * sin_d, sin_i], dtype=float)
        moment_am2 = mag_amp_apm * volume_m3 * moment_unit

        center_xyz_m = lonlatr_to_xyz_m(
            np.array([lat_c_deg]),
            np.array([lon_c_deg]),
            np.array([r_centroid_m]),
        )[0]

        bodies.append(
            {
                "lat_center_deg": lat_c_deg,
                "lon_center_deg": lon_c_deg,
                "volume_m3": volume_m3,
                "depth_top_km": depth_top_km,
                "depth_bot_km": depth_bot_km,
                "center_x_m": center_xyz_m[0],
                "center_y_m": center_xyz_m[1],
                "center_z_m": center_xyz_m[2],
                "mx_am2": moment_am2[0],
                "my_am2": moment_am2[1],
                "mz_am2": moment_am2[2],
            }
        )
        idx += 7
    return bodies


def make_surface_interp(
    lat_axis_deg: np.ndarray,
    lon_axis_deg: np.ndarray,
    elev_km: np.ndarray,
) -> RegularGridInterpolator:
    return RegularGridInterpolator(
        (lat_axis_deg, lon_axis_deg),
        elev_km,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )


def select_seed_points(
    lon_grid_deg: np.ndarray,
    lat_grid_deg: np.ndarray,
    btot_grid_nt: np.ndarray,
    surface_radius_km: np.ndarray,
    east: np.ndarray,
    north: np.ndarray,
    up: np.ndarray,
    origin_xyz_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    xyz_surface_m = lonlatr_to_xyz_m(lat_grid_deg.ravel(), lon_grid_deg.ravel(), surface_radius_km.ravel() * 1000.0)
    x_km, y_km, _z_km = global_to_local_km(xyz_surface_m, origin_xyz_m, east, north, up)
    points_xy_km = np.column_stack([x_km, y_km])
    btot_flat = btot_grid_nt.ravel()
    order = np.argsort(btot_flat)[::-1]

    chosen_xyz: list[np.ndarray] = []
    chosen_btot: list[float] = []
    for idx_flat in order:
        if not np.isfinite(btot_flat[idx_flat]):
            continue
        candidate_xy = points_xy_km[idx_flat]
        if chosen_xyz:
            chosen_xy = np.column_stack(
                global_to_local_km(np.vstack(chosen_xyz), origin_xyz_m, east, north, up)[:2]
            )
            dist = np.linalg.norm(chosen_xy - candidate_xy.reshape(1, 2), axis=1)
            if np.any(dist < MIN_SEED_SEPARATION_KM):
                continue
        lon_deg = lon_grid_deg.ravel()[idx_flat]
        lat_deg = lat_grid_deg.ravel()[idx_flat]
        radius_km = surface_radius_km.ravel()[idx_flat] + SEED_ALTITUDE_KM
        chosen_xyz.append(lonlatr_to_xyz_m(np.array([lat_deg]), np.array([lon_deg]), np.array([radius_km * 1000.0]))[0])
        chosen_btot.append(float(btot_flat[idx_flat]))
        if len(chosen_xyz) >= NUM_SEEDS:
            break

    if not chosen_xyz:
        raise RuntimeError("Could not select any field-line seeds from the surface output.")
    return np.vstack(chosen_xyz), np.asarray(chosen_btot, dtype=float)


def magnetic_field_xyz_nt(
    positions_m: np.ndarray,
    source_xyz_m: np.ndarray,
    source_moment_am2: np.ndarray,
) -> np.ndarray:
    delta = positions_m[:, None, :] - source_xyz_m[None, :, :]
    r2 = np.sum(delta * delta, axis=2)
    r2 = np.maximum(r2, (SOURCE_SOFTENING_KM * 1000.0) ** 2)
    inv_r = 1.0 / np.sqrt(r2)
    inv_r3 = inv_r / r2
    inv_r5 = inv_r3 / r2
    mdotr = np.sum(source_moment_am2[None, :, :] * delta, axis=2)
    field_t = MU0_OVER_4PI * np.sum(
        3.0 * mdotr[:, :, None] * delta * inv_r5[:, :, None]
        - source_moment_am2[None, :, :] * inv_r3[:, :, None],
        axis=1,
    )
    return field_t * 1.0e9


def compute_btot_slice_nt(
    lat_grid_deg: np.ndarray,
    lon_grid_deg: np.ndarray,
    radius_km_grid: np.ndarray,
    source_xyz_m: np.ndarray,
    source_moment_am2: np.ndarray,
) -> np.ndarray:
    xyz_m = lonlatr_to_xyz_m(
        lat_grid_deg.ravel(),
        lon_grid_deg.ravel(),
        radius_km_grid.ravel() * 1000.0,
    )
    field_nt = magnetic_field_xyz_nt(xyz_m, source_xyz_m, source_moment_am2)
    return np.linalg.norm(field_nt, axis=1).reshape(lat_grid_deg.shape)


def extract_contour_segments(
    lon_grid_deg: np.ndarray,
    lat_grid_deg: np.ndarray,
    data_grid: np.ndarray,
    levels_nt: np.ndarray,
) -> list[tuple[float, list[np.ndarray]]]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    cs = ax.contour(lon_grid_deg, lat_grid_deg, data_grid, levels=levels_nt)
    segments: list[tuple[float, list[np.ndarray]]] = []
    for level, segs in zip(cs.levels, cs.allsegs, strict=False):
        keep = [seg.copy() for seg in segs if seg.shape[0] >= 2]
        segments.append((float(level), keep))
    plt.close(fig)
    return segments


def altitude_for_btot_target_km(
    btot_volume_nt: np.ndarray,
    altitude_levels_km: np.ndarray,
    target_nt: float,
) -> np.ndarray:
    nalt, nlat, nlon = btot_volume_nt.shape
    alt_surface = np.full((nlat, nlon), np.nan, dtype=float)
    for iy in range(nlat):
        for ix in range(nlon):
            column = btot_volume_nt[:, iy, ix]
            if not np.all(np.isfinite(column)):
                continue
            for k in range(nalt - 1):
                v0 = column[k] - target_nt
                v1 = column[k + 1] - target_nt
                if v0 == 0.0:
                    alt_surface[iy, ix] = altitude_levels_km[k]
                    break
                if v0 * v1 <= 0.0:
                    a0 = altitude_levels_km[k]
                    a1 = altitude_levels_km[k + 1]
                    if column[k + 1] == column[k]:
                        alt_surface[iy, ix] = a0
                    else:
                        frac = (target_nt - column[k]) / (column[k + 1] - column[k])
                        alt_surface[iy, ix] = a0 + frac * (a1 - a0)
                    break
    return alt_surface


def smooth_iso_altitude_km(iso_altitude_km: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    valid_mask = np.isfinite(iso_altitude_km)
    if not np.any(valid_mask):
        return np.full_like(iso_altitude_km, np.nan), valid_mask

    labels, label_count = ndimage.label(valid_mask)
    keep_mask = np.zeros_like(valid_mask, dtype=bool)
    for label_id in range(1, label_count + 1):
        component_mask = labels == label_id
        if int(component_mask.sum()) >= ISO_MIN_COMPONENT_PIXELS:
            keep_mask |= component_mask
    if not np.any(keep_mask):
        return np.full_like(iso_altitude_km, np.nan), keep_mask

    values = np.where(keep_mask, iso_altitude_km, 0.0)
    weights = keep_mask.astype(float)
    smooth_values = ndimage.gaussian_filter(values, ISO_SMOOTH_SIGMA_PX, mode="nearest")
    smooth_weights = ndimage.gaussian_filter(weights, ISO_SMOOTH_SIGMA_PX, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        smooth_alt = smooth_values / smooth_weights

    smooth_keep_mask = smooth_weights >= ISO_MASK_KEEP_THRESHOLD
    if not np.any(smooth_keep_mask):
        smooth_keep_mask = keep_mask.copy()

    # Taper the surface height to zero at the edge so it merges cleanly into the WAC plane.
    edge_distance_px = ndimage.distance_transform_edt(smooth_keep_mask)
    max_edge_distance_px = float(np.nanmax(edge_distance_px[smooth_keep_mask]))
    taper_width_px = max(
        1.0,
        min(
            max_edge_distance_px,
            max(2.0 * ISO_SMOOTH_SIGMA_PX, ISO_EDGE_TAPER_START * max_edge_distance_px),
        ),
    )
    taper_fraction = np.clip(edge_distance_px / taper_width_px, 0.0, 1.0)
    taper_fraction = taper_fraction * taper_fraction * (3.0 - 2.0 * taper_fraction)

    smooth_alt = np.where(smooth_keep_mask, smooth_alt * taper_fraction, np.nan)
    smooth_alt = np.clip(smooth_alt, 0.0, None)
    return smooth_alt, smooth_keep_mask


def trace_field_lines(
    seeds_xyz_m: np.ndarray,
    source_xyz_m: np.ndarray,
    source_moment_am2: np.ndarray,
    surface_elev_interp: RegularGridInterpolator,
    lat_bounds_deg: tuple[float, float],
    lon_bounds_deg: tuple[float, float],
    lon_ref_deg: np.ndarray,
) -> list[np.ndarray]:
    seed_count = seeds_xyz_m.shape[0]
    directions = np.tile(np.array([-1.0, 1.0], dtype=float), seed_count)
    positions = np.repeat(seeds_xyz_m, 2, axis=0).astype(float)
    arclength_km = np.zeros(positions.shape[0], dtype=float)
    active = np.ones(positions.shape[0], dtype=bool)
    traces = [[positions[i].copy()] for i in range(positions.shape[0])]

    max_arc_km = TRACE_MAX_ARCLENGTH_KM
    source_centroid_m = np.mean(source_xyz_m, axis=0)

    def inside_window(lon_deg: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
        lon_deg = match_lon_domain(lon_deg, lon_ref_deg)
        lon0, lon1 = lon_bounds_deg
        lon_rel = np.array([angle_distance_deg(lon, lon0) for lon in lon_deg], dtype=float)
        lon_span = angle_distance_deg(lon1, lon0)
        return (
            (lat_deg >= lat_bounds_deg[0])
            & (lat_deg <= lat_bounds_deg[1])
            & (lon_rel >= 0.0)
            & (lon_rel <= lon_span)
        )

    while np.any(active):
        active_idx = np.flatnonzero(active)
        xyz_m = positions[active_idx]
        lon_deg, lat_deg, radius_m = xyz_to_lonlatr_m(xyz_m)
        lon_deg = match_lon_domain(lon_deg, lon_ref_deg)
        elev_km = radius_m / 1000.0 - MOON_MEAN_RADIUS_KM
        surf_km = surface_elev_interp(np.column_stack([lat_deg, lon_deg]))

        keep = (
            inside_window(lon_deg, lat_deg)
            & np.isfinite(surf_km)
            & (elev_km >= surf_km - SURFACE_HIT_BUFFER_KM)
            & (elev_km >= TRACE_ALTITUDE_MIN_KM)
            & (elev_km <= TRACE_ALTITUDE_MAX_KM)
            & (arclength_km[active_idx] < max_arc_km)
        )
        active[active_idx[~keep]] = False
        active_idx = np.flatnonzero(active)
        if active_idx.size == 0:
            break

        xyz_m = positions[active_idx]
        dist_to_centroid_km = np.linalg.norm(xyz_m - source_centroid_m.reshape(1, 3), axis=1) / 1000.0
        step_km = np.clip(
            TRACE_STEP_DISTANCE_SCALE * dist_to_centroid_km,
            TRACE_MIN_STEP_KM,
            TRACE_MAX_STEP_KM,
        )
        step_km = np.minimum(step_km, max_arc_km - arclength_km[active_idx])
        step_m = step_km * 1000.0
        sign = directions[active_idx].reshape(-1, 1)

        k1 = magnetic_field_xyz_nt(xyz_m, source_xyz_m, source_moment_am2)
        k1_norm = np.linalg.norm(k1, axis=1)
        valid = np.isfinite(k1_norm) & (k1_norm >= TRACE_MIN_FIELD_NT)
        active[active_idx[~valid]] = False
        active_idx = active_idx[valid]
        if active_idx.size == 0:
            break

        xyz_m = positions[active_idx]
        step_m = step_m[valid]
        sign = sign[valid]
        k1_dir = sign * (k1[valid] / k1_norm[valid][:, None])

        def rk_dir(pos_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            field = magnetic_field_xyz_nt(pos_m, source_xyz_m, source_moment_am2)
            norm = np.linalg.norm(field, axis=1)
            return field, norm

        k2_field, k2_norm = rk_dir(xyz_m + 0.5 * step_m[:, None] * k1_dir)
        valid = np.isfinite(k2_norm) & (k2_norm >= TRACE_MIN_FIELD_NT)
        active[active_idx[~valid]] = False
        active_idx = active_idx[valid]
        if active_idx.size == 0:
            break

        xyz_m = positions[active_idx]
        step_m = step_m[valid]
        sign = sign[valid]
        k1_dir = k1_dir[valid]
        k2_dir = sign * (k2_field[valid] / k2_norm[valid][:, None])

        k3_field, k3_norm = rk_dir(xyz_m + 0.5 * step_m[:, None] * k2_dir)
        valid = np.isfinite(k3_norm) & (k3_norm >= TRACE_MIN_FIELD_NT)
        active[active_idx[~valid]] = False
        active_idx = active_idx[valid]
        if active_idx.size == 0:
            break

        xyz_m = positions[active_idx]
        step_m = step_m[valid]
        sign = sign[valid]
        k1_dir = k1_dir[valid]
        k2_dir = k2_dir[valid]
        k3_dir = sign * (k3_field[valid] / k3_norm[valid][:, None])

        k4_field, k4_norm = rk_dir(xyz_m + step_m[:, None] * k3_dir)
        valid = np.isfinite(k4_norm) & (k4_norm >= TRACE_MIN_FIELD_NT)
        active[active_idx[~valid]] = False
        active_idx = active_idx[valid]
        if active_idx.size == 0:
            break

        xyz_m = positions[active_idx]
        step_m = step_m[valid]
        k1_dir = k1_dir[valid]
        k2_dir = k2_dir[valid]
        k3_dir = k3_dir[valid]
        k4_dir = sign[valid] * (k4_field[valid] / k4_norm[valid][:, None])

        new_xyz_m = xyz_m + (step_m[:, None] / 6.0) * (k1_dir + 2.0 * k2_dir + 2.0 * k3_dir + k4_dir)
        positions[active_idx] = new_xyz_m
        arclength_km[active_idx] += step_m / 1000.0

        for trace_index, point_xyz_m in zip(active_idx, new_xyz_m, strict=False):
            traces[trace_index].append(point_xyz_m.copy())

    field_lines: list[np.ndarray] = []
    for seed_index in range(seed_count):
        backward = np.asarray(traces[2 * seed_index], dtype=float)
        forward = np.asarray(traces[2 * seed_index + 1], dtype=float)
        parts = []
        if backward.shape[0] >= 2:
            parts.append(backward[::-1][:-1])
        if forward.shape[0] >= 2:
            parts.append(forward)
        if parts:
            line = np.vstack(parts)
            if line.shape[0] >= 2:
                field_lines.append(line)
    return field_lines


if not INPUT_IN.is_file():
    raise FileNotFoundError(f"Missing input file: {INPUT_IN}")

if not INPUT_BRTP.is_file():
    raise FileNotFoundError(
        f"Missing BRTP output: {INPUT_BRTP}. Run run_reiner_gamma_table4_surface_direct.sh first."
    )

try:
    import plotly.graph_objects as go
except Exception as exc:  # pragma: no cover - environment-specific dependency
    raise RuntimeError(
        "plotly is required. Run this script with /Users/danywaller/code/venvs/gravmagpy/bin/python."
    ) from exc

bodies = parse_body_models(INPUT_IN)

plot_lon_deg = np.arange(PLOT_LON_MIN_DEG, PLOT_LON_MAX_DEG + 0.5 * PLOT_STEP_DEG, PLOT_STEP_DEG)
plot_lat_deg = np.arange(PLOT_LAT_MIN_DEG, PLOT_LAT_MAX_DEG + 0.5 * PLOT_STEP_DEG, PLOT_STEP_DEG)
lon_grid_deg, lat_grid_deg = np.meshgrid(plot_lon_deg, plot_lat_deg)
wac_norm, wac_mode = load_wac_patch(lat_grid_deg, lon_grid_deg)

source_xyz_m = np.array([[b["center_x_m"], b["center_y_m"], b["center_z_m"]] for b in bodies], dtype=float)
source_moment_am2 = np.array([[b["mx_am2"], b["my_am2"], b["mz_am2"]] for b in bodies], dtype=float)

altitude_levels_km = np.asarray(ALTITUDE_LEVELS_KM, dtype=float)
btot_volume_nt = np.empty((altitude_levels_km.size, lat_grid_deg.shape[0], lon_grid_deg.shape[1]), dtype=float)
for ia, altitude_km in enumerate(altitude_levels_km):
    radius_slice_km = MOON_MEAN_RADIUS_KM + altitude_km + np.zeros_like(lat_grid_deg)
    btot_volume_nt[ia] = compute_btot_slice_nt(
        lat_grid_deg,
        lon_grid_deg,
        radius_slice_km,
        source_xyz_m,
        source_moment_am2,
    )
iso_altitude_km = altitude_for_btot_target_km(btot_volume_nt, altitude_levels_km, ISO_BTOT_NT)
if not np.any(np.isfinite(iso_altitude_km)):
    raise RuntimeError(
        f"No Btot={ISO_BTOT_NT:.0f} nT crossing found between "
        f"{altitude_levels_km[0]:.1f} and {altitude_levels_km[-1]:.1f} km."
    )
iso_altitude_smooth_km, iso_keep_mask = smooth_iso_altitude_km(iso_altitude_km)
iso_absolute_km = iso_altitude_smooth_km
iso_color = np.where(np.isfinite(iso_altitude_smooth_km), iso_altitude_smooth_km, np.nan)
if not np.any(np.isfinite(iso_altitude_smooth_km)):
    raise RuntimeError("Smoothed iso-surface became empty after component cleanup.")
iso_cmin = float(np.nanmin(iso_color))
iso_cmax = float(np.nanmax(iso_color))
if iso_cmax <= iso_cmin:
    iso_cmax = iso_cmin + 1.0e-6

src_lon_deg, src_lat_deg, src_radius_m = xyz_to_lonlatr_m(source_xyz_m)
src_z_km = src_radius_m / 1000.0 - MOON_MEAN_RADIUS_KM
plane_z_km = np.zeros_like(lat_grid_deg)

fig = go.Figure()
fig.add_trace(
    go.Surface(
        x=lon_grid_deg[::SURFACE_DECIMATE, ::SURFACE_DECIMATE],
        y=lat_grid_deg[::SURFACE_DECIMATE, ::SURFACE_DECIMATE],
        z=plane_z_km[::SURFACE_DECIMATE, ::SURFACE_DECIMATE],
        surfacecolor=wac_norm[::SURFACE_DECIMATE, ::SURFACE_DECIMATE],
        colorscale=[[0.0, "rgb(0,0,0)"], [1.0, "rgb(255,255,255)"]],
        cmin=0.0,
        cmax=1.0,
        opacity=SURFACE_OPACITY,
        showscale=False,
        hovertemplate=(
            "Longitude %{x:.3f}°<br>"
            "Latitude %{y:.3f}°<br>"
            "Altitude %{z:.2f} km<br>"
            "WAC plane<extra></extra>"
        ),
        lighting=dict(ambient=0.92, diffuse=0.65, roughness=0.95, specular=0.05),
        lightposition=dict(x=150, y=250, z=500),
        name="WAC plane",
    )
)
fig.add_trace(
    go.Surface(
        x=lon_grid_deg,
        y=lat_grid_deg,
        z=iso_absolute_km,
        surfacecolor=iso_color,
        colorscale="Turbo",
        cmin=iso_cmin,
        cmax=iso_cmax,
        opacity=ISO_SURFACE_OPACITY,
        showscale=True,
        colorbar=dict(title=f"Altitude of {ISO_BTOT_NT:.0f} nT [km]", len=0.65),
        hovertemplate=(
            "Longitude %{x:.3f}°<br>"
            "Latitude %{y:.3f}°<br>"
            "Altitude %{z:.2f} km<br>"
            f"Btot = {ISO_BTOT_NT:.0f} nT<extra></extra>"
        ),
        name=f"{ISO_BTOT_NT:.0f} nT surface",
    )
)

fig.add_trace(
    go.Scatter3d(
        x=src_lon_deg,
        y=src_lat_deg,
        z=np.zeros_like(src_z_km),
        mode="markers",
        marker=dict(size=SOURCE_MARKER_SIZE, color="black"),
        name="Source centers",
        hovertemplate="Source center projected to altitude 0 plane<extra></extra>",
    )
)

z_min = 0.0
z_max = float(np.nanmax(np.where(np.isfinite(iso_absolute_km), iso_absolute_km, 0.0)))

fig.update_layout(
    title="Reiner Gamma equivalent source magnetic field",
    width=PLOTLY_WIDTH_PX,
    height=PLOTLY_HEIGHT_PX,
    margin=dict(l=0, r=0, b=0, t=45),
    legend=dict(x=0.02, y=0.98),
    scene=dict(
        xaxis=dict(
            title="Longitude [°]",
            showgrid=False,
            zeroline=False,
            showbackground=True,
            backgroundcolor="rgba(255,255,255,0.04)",
            showspikes=False,
        ),
        yaxis=dict(
            title="Latitude [°]",
            showgrid=False,
            zeroline=False,
            showbackground=True,
            backgroundcolor="rgba(255,255,255,0.04)",
            showspikes=False,
        ),
        zaxis=dict(
            title="Altitude [km]",
            range=[z_min, z_max],
            showgrid=True,
            zeroline=False,
            showbackground=False,
            gridcolor="rgba(60,60,60,0.55)",
            gridwidth=2,
            showspikes=False,
        ),
        aspectmode="cube",
        camera=dict(eye=dict(x=1.55, y=-1.45, z=0.95)),
    ),
)
OUTPUT_HTML.parent.mkdir(parents=True, exist_ok=True)
fig.write_html(str(OUTPUT_HTML), include_plotlyjs=True)

print(f"Input IN: {INPUT_IN}")
print(f"Input BRTP: {INPUT_BRTP}")
print(f"Wrote {ISO_BTOT_NT:.0f} nT surface HTML: {OUTPUT_HTML}")
print("Surface mode: WAC plane at altitude 0 km")
print(f"WAC mode: {wac_mode}")
print(f"Isosurface target [nT]: {ISO_BTOT_NT}")
print(f"Altitude levels [km]: {ALTITUDE_LEVELS_KM[0]} .. {ALTITUDE_LEVELS_KM[-1]} step {ALTITUDE_LEVELS_KM[1]-ALTITUDE_LEVELS_KM[0]}")
print(f"Finite {ISO_BTOT_NT:.0f} nT cells: {int(np.isfinite(iso_altitude_km).sum())}")
