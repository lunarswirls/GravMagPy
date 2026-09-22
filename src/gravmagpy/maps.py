"""georeferenced context maps alongside equivalent-source solutions"""

from pathlib import Path
import re

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import Affine, from_bounds
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import Window

from .equivalent import equivalent_field_grid, equivalent_source_table
from .geometry import geographic_grid
from .utils.paths import artifact_path


def geographic_crs(radius_km=1737.4):
    """return a planet-centered spherical longitude/latitude crs, not earth wgs84"""
    if not np.isfinite(radius_km) or radius_km <= 0:
        raise ValueError("radius_km must be positive and finite")
    return CRS.from_string(f"+proj=longlat +R={radius_km*1000:.12g} +no_defs")


def read_geotiff(path, *, bounds, shape=(512, 512), radius_km=1737.4, band=1,
                 resampling="bilinear", source_crs=None, nodata=None, scale=None, offset=None):
    """read a bounded, north-up geographic view using the raster's actual georeferencing

    projected/rotated rasters are warped to the requested geographic grid
    nodata/masks are honored; scale/offset default to the selected band metadata
    source_crs is an explicit override for absent or known-incorrect metadata
    categorical maps should use resampling='nearest'
    """
    grid = geographic_grid(bounds, shape)
    destination = geographic_crs(radius_km)
    path = Path(path).expanduser().resolve()
    if resampling not in ("nearest", "bilinear", "cubic", "average"):
        raise ValueError("resampling must be nearest, bilinear, cubic, or average")
    west, south, east, north = grid["bounds"]
    rows, columns = grid["shape"]
    dx = (east-west)/columns
    data = np.full((rows, columns), np.nan)
    with rasterio.open(path) as source:
        crs = CRS.from_user_input(source_crs) if source_crs is not None else source.crs
        if crs is None:
            raise ValueError(f"{path.name} has no crs; supply source_crs explicitly")
        spheroid = re.search(r'(?:SPHEROID|ELLIPSOID)\["[^"]+",\s*([-+0-9.eE]+)', crs.to_wkt())
        if spheroid is None or not np.isclose(float(spheroid.group(1)), radius_km*1000, rtol=0.01):
            raise ValueError("GeoTIFF celestial-body radius differs from the model; check radius_km and source_crs")
        if not isinstance(band, (int, np.integer)) or not 1 <= band <= source.count:
            raise ValueError(f"band must be an integer from 1 to {source.count}")
        factor = source.scales[band-1] if scale is None else float(scale)
        shift = source.offsets[band-1] if offset is None else float(offset)
        if not np.isfinite([factor, shift]).all():
            raise ValueError("scale and offset must be finite")
        # split at both the antimeridian and the source raster's longitude storage seam
        cycles = np.floor((grid["longitude_deg"]+180)/360).astype(int)
        translations = np.zeros(columns, dtype=int)
        if crs.is_geographic:
            center = (source.bounds.left+source.bounds.right)/2
            canonical = grid["longitude_deg"]-360*cycles
            translations = np.floor((canonical-center+180)/360).astype(int)
        breaks = np.flatnonzero((np.diff(cycles) != 0) | (np.diff(translations) != 0))+1
        edges = np.concatenate(([0], breaks, [columns]))
        for start, stop in zip(edges[:-1], edges[1:]):
            cycle = cycles[start]
            left, right = west+start*dx-360*cycle, west+stop*dx-360*cycle
            transform = from_bounds(left, south, right, north, stop-start, rows)
            source_transform = Affine.translation(360*translations[start], 0)*source.transform
            native = transform_bounds(destination, crs, left, south, right, north, densify_pts=41)
            if not np.isfinite(native).all():
                raise ValueError("requested region cannot be transformed into the GeoTIFF projection")
            xmin, ymin, xmax, ymax = native
            corners = [~source_transform * (x, y) for x in (xmin, xmax) for y in (ymin, ymax)]
            pixels = np.asarray(corners)
            col0, row0 = np.maximum(0, np.floor(pixels.min(axis=0)).astype(int)-3)
            col1, row1 = np.minimum([source.width, source.height], np.ceil(pixels.max(axis=0)).astype(int)+3)
            if col1 <= col0 or row1 <= row0:
                continue
            if (col1-col0)*(row1-row0) > 16_000_000:
                raise ValueError("source window exceeds sixteen million pixels; use a smaller region or downsampled GeoTIFF")
            window = Window(int(col0), int(row0), int(col1-col0), int(row1-row0))
            # read explicit validity masks before warping; gdal's virtual warp may discard them
            values = source.read(band, window=window, masked=True, out_dtype="float64").filled(np.nan)
            for invalid in (source.nodatavals[band-1], nodata):
                if invalid is not None:
                    values[values == invalid] = np.nan
            values[~np.isfinite(values)] = np.nan
            warped = np.full((rows, stop-start), np.nan)
            reproject(values, warped, src_transform=source_transform*Affine.translation(col0, row0),
                      src_crs=crs, src_nodata=np.nan, dst_transform=transform, dst_crs=destination,
                      dst_nodata=np.nan, resampling=Resampling[resampling], warp_mem_limit=64)
            data[:, start:stop] = warped*factor + shift
        source_wkt = crs.to_wkt()
        unit = source.units[band-1]
    data[~np.isfinite(data)] = np.nan
    if not np.isfinite(data).any():
        raise ValueError(f"{path.name} has no valid pixels in the requested region")
    return {**grid, "data": np.ma.masked_invalid(data), "path": path, "band": int(band),
            "crs": destination, "source_crs": source_wkt, "scale": float(factor), "offset": float(shift),
            "unit": unit, "transform": from_bounds(west, south, east, north, columns, rows)}


def _color_limits(values, vmin=None, vmax=None, symmetric=False):
    values = np.ma.masked_invalid(values).compressed()
    if not len(values):
        raise ValueError("cannot plot a panel with no finite values")
    low, high = np.percentile(values, [2, 98])
    if symmetric:
        high = max(abs(low), abs(high))
        low = -high
    low = float(low if vmin is None else vmin)
    high = float(high if vmax is None else vmax)
    if not np.isfinite([low, high]).all() or low > high:
        raise ValueError("color limits must be finite and increasing")
    if low == high:
        delta = max(abs(low)*0.01, 1e-12)
        low, high = low-delta, high+delta
    return low, high


def plot_equivalent_maps(solution, rasters, image_path=None, *, bounds, altitude_km=30.0,
                         shape=(241, 241), radius_km=None, components=("btot",), source_layer=1,
                         show_source_panel=True, show_sources=True, contours=True, contour_levels=6,
                         ncols=2, title=None, dpi=200):
    """plot equivalent fields alongside georeferenced context panels

    solution is a fit result, source dataframe, or dipoles.csv path
    rasters is a list of dictionaries with path, title, cmap, unit, vmin/vmax,
    plus optional read_geotiff band/resampling/source_crs/nodata/scale/offset
    all source layers contribute to predicted fields; source_layer selects only
    which layer's moments and surface locations are drawn
    """
    from .plotting import plt
    from matplotlib.ticker import FuncFormatter

    if isinstance(components, str):
        components = (components,)
    if not components or len(set(components)) != len(components) or set(components)-{"bx", "by", "bz", "btot"}:
        raise ValueError("components must select distinct bx, by, bz, and/or btot")
    if not isinstance(ncols, (int, np.integer)) or not 1 <= ncols <= 4:
        raise ValueError("ncols must be an integer from 1 to 4")
    if not isinstance(contour_levels, (int, np.integer)) or not 2 <= contour_levels <= 20:
        raise ValueError("contour_levels must be an integer from 2 to 20")
    if image_path is None:
        source_path = None
        if isinstance(solution, (str, Path)):
            source_path = solution
        elif isinstance(solution, dict):
            source_path = solution.get("source_path")
        if source_path is None:
            raise ValueError("image_path is required for a solution without a source csv path")
        image_path = artifact_path(source_path, kind="figs", suffix=f"_context_{altitude_km:g}km.png")
    image_path = Path(image_path).expanduser().resolve()
    for spec in rasters:
        if Path(spec["path"]).expanduser().resolve() == image_path:
            raise ValueError("image output must not overwrite an input GeoTIFF")
    if isinstance(solution, (str, Path)) and Path(solution).expanduser().resolve() == image_path:
        raise ValueError("image output must not overwrite the source csv")
    grid = equivalent_field_grid(solution, bounds=bounds, altitude_km=altitude_km, shape=shape, radius_km=radius_km)
    radius_km = grid["radius_km"]
    sources = equivalent_source_table(solution)
    if "layer_id" in sources:
        sources = sources.loc[sources.layer_id == source_layer]
    elif source_layer != 1:
        raise ValueError("source_layer must be 1 for a table without layer_id")
    if not len(sources):
        raise ValueError(f"source layer {source_layer} contains no sources")
    west, south, east, north = grid["bounds"]
    center = (west+east)/2
    source_lon = center + (sources.lon_deg.to_numpy()-center+180) % 360 - 180
    source_lat = sources.lat_deg.to_numpy()
    inside = (source_lon >= west) & (source_lon <= east) & (source_lat >= south) & (source_lat <= north)
    moment = np.linalg.norm(sources[["mx_am2", "my_am2", "mz_am2"]].to_numpy(), axis=1)
    panels = []
    labels = {"bx": "Bx", "by": "By", "bz": "Bz", "btot": "Btot"}
    for component in components:
        values = grid[f"{component}_nt"]
        low, high = _color_limits(values, vmin=0 if component == "btot" else None, symmetric=component != "btot")
        panels.append({"kind": "field", "data": values, "title": f"Equivalent {labels[component]} | {altitude_km:g} km",
                       "unit": "nT", "cmap": "viridis" if component == "btot" else "RdBu_r", "vmin": low, "vmax": high})
    if show_source_panel:
        if not inside.any():
            raise ValueError("no source locations in the map bounds; disable show_source_panel or expand bounds")
        low, high = _color_limits(moment[inside], vmin=0)
        depth = radius_km-sources.radius_m.to_numpy()/1000
        depth_label = f"{depth.min():g}" if np.isclose(depth.min(), depth.max()) else f"{depth.min():g}–{depth.max():g}"
        panels.append({"kind": "sources", "title": f"Source moments | layer {source_layer}, {depth_label} km depth",
                       "unit": "moment magnitude [A m²]", "cmap": "magma", "vmin": low, "vmax": high})
    read_keys = {"band", "resampling", "source_crs", "nodata", "scale", "offset"}
    for spec in rasters:
        unknown = set(spec) - read_keys - {"path", "title", "cmap", "unit", "vmin", "vmax"}
        if unknown:
            raise ValueError(f"unknown raster settings: {sorted(unknown)}")
        raster = read_geotiff(spec["path"], bounds=grid["bounds"], shape=grid["shape"], radius_km=radius_km,
                              **{key: spec[key] for key in read_keys if key in spec})
        low, high = _color_limits(raster["data"], spec.get("vmin"), spec.get("vmax"))
        panels.append({"kind": "raster", "data": raster["data"], "title": spec.get("title", raster["path"].stem),
                       "unit": spec.get("unit", raster["unit"] or "map value"), "cmap": spec.get("cmap", "gray"),
                       "vmin": low, "vmax": high})
    nrows = int(np.ceil(len(panels)/ncols))
    figure, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5.2*nrows), squeeze=False, constrained_layout=True)
    levels = np.linspace(float(grid["btot_nt"].min()), float(grid["btot_nt"].max()), contour_levels+2)[1:-1]
    try:
        for index, (axis, panel) in enumerate(zip(axes.ravel(), panels)):
            axis.set_facecolor("0.92")
            if panel["kind"] == "sources":
                artist = axis.scatter(source_lon[inside], source_lat[inside], c=moment[inside], s=65,
                                      edgecolors="black", linewidths=0.4, cmap=panel["cmap"],
                                      vmin=panel["vmin"], vmax=panel["vmax"])
            else:
                artist = axis.imshow(panel["data"], extent=(west, east, south, north), origin="upper",
                                     interpolation="nearest", cmap=panel["cmap"], vmin=panel["vmin"], vmax=panel["vmax"])
                if show_sources and inside.any():
                    axis.scatter(source_lon[inside], source_lat[inside], marker="o", s=14,
                                 facecolors="none", edgecolors="cyan", linewidths=0.7)
            if contours and panel["kind"] == "raster" and np.ptp(grid["btot_nt"]) > 0:
                contour = axis.contour(grid["longitude_deg"], grid["latitude_deg"], grid["btot_nt"],
                                       levels=levels, colors="white", linewidths=0.7)
                axis.clabel(contour, inline=True, fmt="%.2g", fontsize=7)
            axis.set_xlim(west, east)
            axis.set_ylim(south, north)
            axis.set_aspect(1/max(np.cos(np.deg2rad((south+north)/2)), 0.05))
            axis.xaxis.set_major_formatter(FuncFormatter(lambda value, position: f"{(value+180)%360-180:g}"))
            axis.set_xlabel("longitude [°E]")
            axis.set_ylabel("latitude [°N]")
            axis.set_title(f"({chr(97+index)}) {panel['title']}", fontsize=11)
            figure.colorbar(artist, ax=axis, shrink=0.85, pad=0.025).set_label(panel["unit"])
        for axis in axes.ravel()[len(panels):]:
            axis.set_visible(False)
        notes = []
        if show_sources:
            notes.append(f"cyan circles: source layer {source_layer} projected onto the surface")
        if contours and rasters:
            notes.append(f"white contours: equivalent Btot at {altitude_km:g} km [nT]")
        heading = title or "Equivalent-source solution and georeferenced context maps"
        figure.suptitle(heading + ("\n" + "; ".join(notes) if notes else ""), fontsize=12)
        image_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(image_path, dpi=dpi)
    finally:
        plt.close(figure)
    return image_path
