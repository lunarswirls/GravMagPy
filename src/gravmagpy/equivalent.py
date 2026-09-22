"""python interface to the gravmag sphere equivalent dipole-grid fitter"""

import ctypes
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .forward import load_dipole_backend
from .geometry import geographic_grid, spherical_to_cartesian
from .utils import build_fortran, run_fortran


def equivalent_source_table(solution):
    """read fitted dipoles from a fit result, dataframe, or saved dipoles.csv"""
    if isinstance(solution, dict):
        solution = solution["source_table"]
    table = pd.read_csv(solution) if isinstance(solution, (str, Path)) else pd.DataFrame(solution).copy()
    required = ["lat_deg", "lon_deg", "radius_m", "mx_am2", "my_am2", "mz_am2"]
    if not set(required).issubset(table.columns) or not len(table):
        raise ValueError(f"equivalent sources require a nonempty table with {required}")
    values = table[required].to_numpy(dtype=float)
    if not np.isfinite(values).all() or np.any(np.abs(values[:, 0]) > 90) or np.any(values[:, 2] <= 0):
        raise ValueError("source positions and moments must be finite with valid latitude and positive radius")
    table[required] = values
    return table


def predict_equivalent_field(xyz_km, solution, *, radius_km=None, library=None):
    """evaluate equivalent dipole moments with fortran at new positions in km

    return planet-fixed bx/by/bz in nt; moments are a m2, not magnetization a/m
    observations may lie on or above the reference sphere; csv inputs default
    to the lunar radius unless radius_km is supplied explicitly
    """
    if radius_km is None:
        radius_km = solution.get("radius_km", 1737.4) if isinstance(solution, dict) else 1737.4
    if not np.isfinite(radius_km) or radius_km <= 0:
        raise ValueError("radius_km must be finite and positive")
    table = equivalent_source_table(solution)
    if np.any(table.radius_m.to_numpy() > radius_km*1000 + 1e-6):
        raise ValueError("equivalent sources must be on or below the specified reference surface")
    xyz = np.ascontiguousarray(xyz_km, dtype=np.float64)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError("xyz_km must be a finite nonempty (n, 3) array")
    if np.any(np.linalg.norm(xyz, axis=1) < radius_km-1e-9):
        raise ValueError("map observations must be on or above the reference sphere")
    if max(len(xyz), len(table)) > np.iinfo(np.int32).max:
        raise ValueError("array dimensions exceed the fortran integer range")
    positions = np.ascontiguousarray(spherical_to_cartesian(table.lat_deg, table.lon_deg, table.radius_m/1000))
    moments = np.ascontiguousarray(table[["mx_am2", "my_am2", "mz_am2"]].to_numpy())
    backend = load_dipole_backend(library or build_fortran("orbital"))
    field = np.zeros_like(xyz)
    status = ctypes.c_int()
    backend.dipole_field(len(xyz), len(table), xyz, positions, moments, field, ctypes.byref(status))
    if status.value or not np.isfinite(field).all():
        raise RuntimeError(f"Fortran equivalent field failed with status {status.value}; check source/observer separation")
    return field


def equivalent_field_grid(solution, *, bounds, altitude_km=30.0, shape=(181, 181), radius_km=None, library=None):
    """predict a north-up, constant-altitude map from fitted equivalent sources"""
    if radius_km is None:
        radius_km = solution.get("radius_km", 1737.4) if isinstance(solution, dict) else 1737.4
    if not np.isfinite(altitude_km) or altitude_km < 0:
        raise ValueError("altitude_km must be finite and nonnegative")
    grid = geographic_grid(bounds, shape)
    longitude, latitude = np.meshgrid(grid["longitude_deg"], grid["latitude_deg"])
    xyz = spherical_to_cartesian(latitude.ravel(), longitude.ravel(), radius_km+altitude_km)
    field = predict_equivalent_field(xyz, solution, radius_km=radius_km, library=library).reshape(*grid["shape"], 3)
    grid.update(field=field, btot_nt=np.linalg.norm(field, axis=-1), altitude_km=float(altitude_km), radius_km=float(radius_km))
    grid.update({f"b{axis}_nt": field[..., index] for index, axis in enumerate(("x", "y", "z"))})
    return grid


def fit_equivalent_sources(observations, *, depth_km=10.0, spacing_deg=(1.0, 1.0),
                           layers=1, layer_spacing_km=5.0, padding_deg=(0.0, 0.0),
                           regularization=0.0, max_memory_mib=768.0, output_dir=None, **build_options):
    """fit fixed-grid vector dipole moments using the existing fortran solver

    observations uses the same xyz_km/field_nt/radius_km dictionary as the orbital
    interface; depth_km is below the surface, not below the lowest observation
    spacing_deg and padding_deg are latitude/longitude pairs
    regularization is the native fortran ridge parameter, defaulting here to zero
    additional damping; the solver still includes its fixed diagonal floor
    this solver is unweighted and has no fitted background offsets; nonuniform
    uncertainties are rejected rather than silently discarded
    """
    radius = float(observations["radius_km"])
    xyz = np.asarray(observations["xyz_km"], dtype=float)
    observed = np.asarray(observations["field_nt"], dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 1 or observed.shape != xyz.shape:
        raise ValueError("xyz_km and field_nt must be matching nonempty (n, 3) arrays")
    if not np.isfinite(xyz).all() or not np.isfinite(observed).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError("positions/fields must be finite and radius_km positive")
    if observations.get("frame", "planet_fixed") not in ("SEL", "planet_fixed"):
        raise ValueError("equivalent sources require planet-fixed positions and field vectors")
    sigma = np.broadcast_to(np.asarray(observations.get("sigma_nt", 1.0), dtype=float), observed.shape)
    if not np.isfinite(sigma).all() or np.any(sigma <= 0) or not np.all(sigma == sigma.flat[0]):
        raise ValueError("the fortran dipole-grid fitter requires uniform positive uncertainties")
    radii = np.linalg.norm(xyz, axis=1)
    if np.any(radii <= radius):
        raise ValueError("observations must be above the reference surface")
    scalars = np.asarray([depth_km, layers, layer_spacing_km, regularization, max_memory_mib], dtype=float)
    if not np.isfinite(scalars).all() or depth_km < 0 or layers < 1 or int(layers) != layers or regularization < 0 or max_memory_mib <= 0:
        raise ValueError("invalid depth, layer count, regularization, or memory limit")
    if layer_spacing_km < 0 or (layers > 1 and layer_spacing_km <= 0) or depth_km + (layers-1)*layer_spacing_km >= radius:
        raise ValueError("radial layers must have positive spacing and stay above the planet center")
    spacing, padding = np.asarray(spacing_deg, dtype=float), np.asarray(padding_deg, dtype=float)
    if spacing.shape != (2,) or padding.shape != (2,) or not np.isfinite([spacing, padding]).all() or np.any(spacing <= 0) or np.any(padding < 0):
        raise ValueError("spacing_deg must be a positive pair and padding_deg a nonnegative pair")
    latitude = np.rad2deg(np.arctan2(xyz[:, 2], np.hypot(xyz[:, 0], xyz[:, 1])))
    longitude = np.rad2deg(np.arctan2(xyz[:, 1], xyz[:, 0]))
    # mirror the fitter's grid sizing before its source arrays are allocated
    angle = np.deg2rad(longitude)
    reference = np.rad2deg(np.arctan2(np.sin(angle).sum(), np.cos(angle).sum()))
    unwrapped = reference + (longitude-reference+180) % 360 - 180
    low_lat = max(-89.75, latitude.min()-padding[0])
    high_lat = min(89.75, latitude.max()+padding[0])
    nlat = max(1, int(np.ceil((high_lat-low_lat)/spacing[0]))+1)
    nlon = max(1, int(np.ceil((np.ptp(unwrapped)+2*padding[1])/spacing[1]))+1)
    if low_lat + (nlat-1)*spacing[0] >= 90 or (nlon-1)*spacing[1] >= 360:
        raise ValueError("equivalent grid crosses a pole or repeats longitude; reduce padding/spacing or select a local region")
    nsource = nlat*nlon*int(layers)
    nrow, ncolumn = 3*len(xyz), 3*nsource
    memory_mib = (nrow*ncolumn + ncolumn*ncolumn + 4*ncolumn + 2*nrow)*8/1024**2
    if max(nrow, ncolumn, nsource) > np.iinfo(np.int32).max or memory_mib > max_memory_mib:
        raise ValueError(f"equivalent fit needs approximately {memory_mib:.3f} MiB; increase spacing or the memory limit")
    below_observations = depth_km + radii.min() - radius
    with tempfile.TemporaryDirectory(prefix="gravmagpy-equivalent-") as temporary:
        temporary = Path(temporary)
        # short internal paths avoid the fixed-length fortran cli buffers and csv path-list limit
        pd.DataFrame({"latitude_deg": latitude, "longitude_deg": longitude, "radius_km": radii,
                      "bx_nt": observed[:, 0], "by_nt": observed[:, 1], "bz_nt": observed[:, 2]}).to_csv(
                          temporary / "observations.csv", index=False)
        run = run_fortran("dipole_grid", [radius, "observations.csv", "predictions.csv", "dipoles.csv",
                          below_observations, *spacing, regularization, int(layers), layer_spacing_km, *padding,
                          max_memory_mib], cwd=temporary, **build_options)
        predictions = pd.read_csv(temporary / "predictions.csv")
        dipoles = pd.read_csv(temporary / "dipoles.csv")
    if len(predictions) != len(xyz) or len(dipoles) != nsource:
        raise RuntimeError("unexpected fortran equivalent-source or observation count")
    predicted = predictions[["bx_pred_nt", "by_pred_nt", "bz_pred_nt"]].to_numpy()
    moments = dipoles[["mx_am2", "my_am2", "mz_am2"]].to_numpy()
    if not np.isfinite(predicted).all() or not np.isfinite(dipoles.to_numpy()).all():
        raise RuntimeError("nonfinite values in equivalent-source output")
    dipoles["depth_km"] = radius - dipoles["radius_m"] / 1000
    if "table" in observations:
        metadata = observations["table"].reset_index(drop=True)
        for key in ("source_file", "source_row", "utc"):
            if key in metadata:
                predictions[key] = metadata[key]
    residual = predicted-observed
    result = {"source_table": dipoles, "predictions": predictions, "predicted_nt": predicted,
              "success": True, "message": "fortran linear solve completed",
              "residual_nt": residual, "moments_a_m2": moments, "radius_km": radius,
              "rmse_nt": float(np.sqrt(np.mean(residual**2))),
              "component_rmse_nt": np.sqrt(np.mean(residual**2, axis=0)).tolist(),
              "grid_shape": (int(layers), nlat, nlon), "depth_km": depth_km,
              "depth_below_min_observation_km": float(below_observations), "regularization": regularization,
              "estimated_memory_mib": memory_mib, "stdout": run.stdout, "stderr": run.stderr}
    if output_dir is not None:
        directory = Path(output_dir).expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        dipoles.to_csv(directory / "dipoles.csv", index=False)
        predictions.to_csv(directory / "predictions.csv", index=False)
        result["output_dir"] = directory
        result["source_path"] = directory / "dipoles.csv"
    return result
