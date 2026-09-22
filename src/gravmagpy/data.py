"""read wake-sorted lunar prospector csv files in the body-fixed sel frame"""

from pathlib import Path

import numpy as np
import pandas as pd

from .geometry import spherical_to_cartesian


def load_lpmag_csv(paths, *, radius_km=1737.4, columns=None, wake_only=True,
                   sigma_nt=1.0, sigma_columns=None, latitude_range=None, longitude_range=None,
                   altitude_range_km=None, max_brms_nt=None):
    """load one or more csv files, preserving row order and per-sample altitude

    default columns: X_SEL/Y_SEL/Z_SEL [km], Bx_SEL/By_SEL/Bz_SEL [nt]
    columns maps lowercase canonical keys to actual headers; use lat_deg,
    lon_deg, altitude_km instead of x_km/y_km/z_km for spherical input
    wake_only uses finite positive wake_t_km and finite wake_rperp_km when
    present, otherwise requires a filename containing 'wake' (preselected data)
    Brms is field variability and is only used as a quality filter, never
    silently interpreted as independent component measurement uncertainties
    """
    if not np.isfinite(radius_km) or radius_km <= 0:
        raise ValueError("radius_km must be finite and positive")
    if isinstance(paths, (str, Path)):
        paths = [paths]
    paths = [Path(path).expanduser().resolve() for path in paths]
    if not paths or len(set(paths)) != len(paths):
        raise ValueError("provide at least one csv path without duplicate paths")
    mapping = {
        "x_km": "X_SEL", "y_km": "Y_SEL", "z_km": "Z_SEL",
        "bx_nt": "Bx_SEL", "by_nt": "By_SEL", "bz_nt": "Bz_SEL", "utc": "utc",
        "brms_nt": "Brms", "wake_t_km": "wake_t_km", "wake_rperp_km": "wake_rperp_km",
    }
    unknown = set(columns or {}) - set(mapping) - {"lat_deg", "lon_deg", "altitude_km"}
    if unknown:
        raise ValueError(f"unknown canonical column keys: {sorted(unknown)}")
    mapping.update(columns or {})
    spherical = all(key in mapping for key in ("lat_deg", "lon_deg", "altitude_km"))
    position_keys = ("lat_deg", "lon_deg", "altitude_km") if spherical else ("x_km", "y_km", "z_km")
    field_keys = ("bx_nt", "by_nt", "bz_nt")
    frames, reports = [], []
    for path in paths:
        frame = pd.read_csv(path)
        required = [mapping[key] for key in (*position_keys, *field_keys)]
        if sigma_columns is not None:
            if len(sigma_columns) != 3:
                raise ValueError("sigma_columns must contain three component uncertainty column names")
            required += list(sigma_columns)
        missing = set(required) - set(frame.columns)
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
        frame = frame.copy()
        frame[required] = frame[required].apply(pd.to_numeric, errors="coerce")
        mask = np.isfinite(frame[required].to_numpy(dtype=float)).all(axis=1)
        if wake_only:
            wake_columns = [mapping[key] for key in ("wake_t_km", "wake_rperp_km")]
            if all(name in frame for name in wake_columns):
                wake = frame[wake_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
                mask &= np.isfinite(wake).all(axis=1) & (wake[:, 0] > 0) & (wake[:, 1] >= 0)
            elif "wake" not in path.stem.lower():
                raise ValueError(f"{path.name}: no wake metadata; use a preselected '*wake*.csv' or wake_only=False")
        if max_brms_nt is not None:
            if mapping["brms_nt"] not in frame:
                raise ValueError(f"{path.name}: Brms column required for max_brms_nt")
            brms = pd.to_numeric(frame[mapping["brms_nt"]], errors="coerce").to_numpy()
            mask &= np.isfinite(brms) & (brms >= 0) & (brms <= max_brms_nt)
        input_rows = len(frame)
        row_ids = np.flatnonzero(mask)
        frame = frame.loc[mask].copy()
        position = frame[[mapping[key] for key in position_keys]].to_numpy(dtype=float)
        if spherical:
            if np.any(np.abs(position[:, 0]) > 90):
                raise ValueError(f"{path.name}: latitude outside [-90, 90]")
            if np.any(position[:, 2] <= 0):
                raise ValueError(f"{path.name}: altitude_km must be positive for orbital samples")
            xyz = spherical_to_cartesian(position[:, 0], position[:, 1], radius_km + position[:, 2])
        else:
            xyz = position
        radii = np.linalg.norm(xyz, axis=1)
        if np.any(radii <= radius_km):
            raise ValueError(f"{path.name}: positions at/below the reference sphere; check km units and radius_km")
        lon = np.rad2deg(np.arctan2(xyz[:, 1], xyz[:, 0]))
        lat = np.rad2deg(np.arctan2(xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1])))
        altitude = radii - radius_km
        sigma = frame[list(sigma_columns)].to_numpy(dtype=float) if sigma_columns else np.broadcast_to(sigma_nt, (len(frame), 3)).copy()
        if not np.isfinite(sigma).all() or np.any(sigma <= 0):
            raise ValueError("component uncertainties must be finite and positive")
        subset = np.ones(len(frame), dtype=bool)
        for interval, values in ((latitude_range, lat), (altitude_range_km, altitude)):
            if interval is not None:
                low, high = interval
                if not np.isfinite([low, high]).all() or low > high:
                    raise ValueError("selection ranges must be finite (lower, upper) pairs")
                subset &= (values >= low) & (values <= high)
        if longitude_range is not None:
            low, high = longitude_range
            if not np.isfinite([low, high]).all() or not (-180 <= low <= 180 and -180 <= high <= 180):
                raise ValueError("longitude_range must use degrees within [-180, 180]")
            subset &= ((lon >= low) & (lon <= high)) if low <= high else ((lon >= low) | (lon <= high))
        normalized = pd.DataFrame({
            "source_file": str(path), "source_row": row_ids,
            "utc": frame[mapping["utc"]].to_numpy() if mapping["utc"] in frame else "",
            "lon_deg": lon, "lat_deg": lat, "altitude_km": altitude,
        })
        for index, key in enumerate(("x_km", "y_km", "z_km")):
            normalized[key] = xyz[:, index]
        for index, key in enumerate(field_keys):
            normalized[key] = frame[mapping[key]].to_numpy()
            normalized[f"sigma_{key}"] = sigma[:, index]
        frames.append(normalized.loc[subset])
        reports.append({"path": str(path), "input_rows": input_rows, "retained_rows": int(subset.sum())})
    table = pd.concat(frames, ignore_index=True)
    if table.empty:
        raise ValueError("no valid observations remain after selection")
    return {
        "xyz_km": table[["x_km", "y_km", "z_km"]].to_numpy(dtype=float),
        "field_nt": table[list(field_keys)].to_numpy(dtype=float),
        "sigma_nt": table[[f"sigma_{key}" for key in field_keys]].to_numpy(dtype=float),
        "radius_km": float(radius_km), "frame": "SEL", "table": table, "report": reports,
    }
