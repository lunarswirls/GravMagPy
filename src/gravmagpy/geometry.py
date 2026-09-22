"""spherical source blocks in a planet-fixed cartesian magnetic frame"""

import numpy as np

parameter_names = (
    "lat_deg", "lon_deg", "lat_width_deg", "lon_width_deg", "depth_top_km", "thickness_km",
    "mx_a_m", "my_a_m", "mz_a_m",
)


def spherical_to_cartesian(lat_deg, lon_deg, radius_km):
    """return (..., 3) planet-fixed xyz positions in km"""
    lat, lon, radius = np.broadcast_arrays(np.deg2rad(lat_deg), np.deg2rad(lon_deg), radius_km)
    return np.stack((radius*np.cos(lat)*np.cos(lon), radius*np.cos(lat)*np.sin(lon), radius*np.sin(lat)), axis=-1)


def geographic_grid(bounds, shape=(181, 181)):
    """return north-up pixel-center axes for (west, south, east, north) degree bounds

    shape is (rows, columns); east < west declares an antimeridian crossing
    longitude remains continuous across that seam, for example 170 to 190
    """
    bounds = np.asarray(bounds, dtype=float)
    dimensions = np.asarray(shape, dtype=float)
    if bounds.shape != (4,) or not np.isfinite(bounds).all():
        raise ValueError("bounds must be finite (west, south, east, north) degrees")
    if dimensions.shape != (2,) or not np.isfinite(dimensions).all() or np.any(dimensions < 2) or np.any(dimensions != np.floor(dimensions)):
        raise ValueError("shape must contain two integer dimensions of at least two")
    if np.prod(dimensions) > 4_000_000:
        raise ValueError("map grids are limited to four million pixels; use a coarser shape")
    west, south, east, north = bounds
    if east < west:
        east += 360
    if not 0 < east-west <= 180 or not -90 <= south < north <= 90:
        raise ValueError("bounds must have positive area, valid latitudes, and at most 180 degrees longitude span")
    rows, columns = map(int, dimensions)
    dx, dy = (east-west)/columns, (north-south)/rows
    return {"bounds": (float(west), float(south), float(east), float(north)), "shape": (rows, columns),
            "longitude_deg": west + (np.arange(columns)+0.5)*dx,
            "latitude_deg": north - (np.arange(rows)+0.5)*dy}


def source_array(sources, radius_km):
    """validate source dictionaries and pack them for the fortran interface"""
    if not np.isfinite(radius_km) or radius_km <= 0:
        raise ValueError("radius_km must be positive and finite")
    if not sources:
        raise ValueError("at least one source is required")
    rows = np.asarray([[source[key] for key in parameter_names] for source in sources], dtype=np.float64)
    if rows.ndim != 2 or rows.shape[1] != 9 or not np.isfinite(rows).all():
        raise ValueError("each source parameter must be a finite scalar")
    if np.any(rows[:, 2:4] <= 0) or np.any(rows[:, 3] > 360):
        raise ValueError("source angular widths must be positive; longitude width cannot exceed 360 degrees")
    if np.any(np.abs(rows[:, 0]) + rows[:, 2] / 2 >= 90):
        raise ValueError("source latitude bounds must stay strictly between the poles")
    if np.any(rows[:, 4] < 0) or np.any(rows[:, 5] <= 0) or np.any(rows[:, 4] + rows[:, 5] >= radius_km):
        raise ValueError("sources need nonnegative top depth, positive thickness, and bottom depth below the radius")
    return np.ascontiguousarray(rows)


def block_volume(sources, radius_km=1737.4):
    """return exact spherical block volumes in cubic meters"""
    rows = source_array(sources, radius_km)
    lat_low = np.deg2rad(rows[:, 0] - rows[:, 2] / 2)
    lat_high = np.deg2rad(rows[:, 0] + rows[:, 2] / 2)
    top = (radius_km - rows[:, 4]) * 1000
    bottom = top - rows[:, 5] * 1000
    return np.deg2rad(rows[:, 3]) * (np.sin(lat_high) - np.sin(lat_low)) * (top**3 - bottom**3) / 3
