"""bounded orbital source geometry fitting with scipy and a fortran forward model"""

from copy import deepcopy
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .forward import predict_field
from .geometry import parameter_names, source_array
from .utils import artifact_path, build_fortran, observation_input_path


def fit_sources(observations, initial_sources, bounds, *, quadrature_order=6, fit_background=False,
                priors=None, loss="linear", f_scale=1.0, max_nfev=200, library=None):
    """fit source geometry/magnetization to all orbital altitudes simultaneously

    bounds is one dictionary per source: {parameter: (lower, upper)}
    omitted parameters stay fixed; the number of bodies is specified by the user
    optional priors use {parameter: (mean, standard_deviation)} per source
    fit_background adds a constant xyz field for each source_file (or each group
    in observations['groups']); these offsets are jointly optimized
    uncertainties are observations['sigma_nt'], in nt, for each vector component
    """
    radius = float(observations["radius_km"])
    source_array(initial_sources, radius)
    xyz = np.asarray(observations["xyz_km"], dtype=float)
    observed = np.asarray(observations["field_nt"], dtype=float)
    if observed.shape != xyz.shape or observed.ndim != 2 or observed.shape[1] != 3 or not np.isfinite(observed).all():
        raise ValueError("field_nt and xyz_km must be finite matching (n, 3) arrays")
    sigma = np.broadcast_to(observations.get("sigma_nt", 1.0), observed.shape)
    if not np.isfinite(sigma).all() or np.any(sigma <= 0):
        raise ValueError("sigma_nt must be finite and positive")
    if len(bounds) != len(initial_sources):
        raise ValueError("provide one bounds dictionary for every source")
    priors = priors if priors is not None else [{} for source in initial_sources]
    if len(priors) != len(initial_sources):
        raise ValueError("provide one priors dictionary for every source")
    indices, initial, lower, upper, names = [], [], [], [], []
    for index, (source, limits) in enumerate(zip(initial_sources, bounds)):
        unknown = set(limits) - set(parameter_names)
        if unknown:
            raise ValueError(f"unknown source parameters: {sorted(unknown)}")
        for key, (low, high) in limits.items():
            value = float(source[key])
            if not np.isfinite([low, high]).all() or not low < high or not low <= value <= high:
                raise ValueError(f"source {index} {key}: need finite lower < upper with initial value inside bounds")
            indices.append((index, key))
            names.append(f"source_{index}.{key}")
            initial.append(value)
            lower.append(low)
            upper.append(high)
        # ensure the entire rectangular bound region describes valid subsurface blocks
        low_source, high_source = dict(source), dict(source)
        for key, (low, high) in limits.items():
            low_source[key], high_source[key] = low, high
        source_array([low_source, high_source], radius)
        max_latitude = max(abs(low_source["lat_deg"]), abs(high_source["lat_deg"]))
        if max_latitude + high_source["lat_width_deg"] / 2 >= 90:
            raise ValueError("latitude and width bounds can cross a pole")
        for key, (mean, std) in priors[index].items():
            if key not in limits or not np.isfinite([mean, std]).all() or std <= 0:
                raise ValueError("priors require fitted parameters with finite means and positive standard deviations")
    if not indices:
        raise ValueError("select at least one source parameter to fit using bounds")
    source_parameter_count = len(initial)
    group_names, group_index = np.array([], dtype=str), None
    background_only_rmse = None
    if fit_background:
        if "groups" in observations:
            groups = np.asarray(observations["groups"]).astype(str)
        elif "table" in observations:
            groups = observations["table"]["source_file"].to_numpy().astype(str)
        else:
            groups = np.full(len(xyz), "all")
        if groups.shape != (len(xyz),):
            raise ValueError("background groups must contain one label per observation")
        group_names, group_index = np.unique(groups, return_inverse=True)
        group_weight = np.zeros((len(group_names), 3))
        group_field = np.zeros_like(group_weight)
        np.add.at(group_weight, group_index, 1 / sigma**2)
        np.add.at(group_field, group_index, observed / sigma**2)
        background_only = (group_field / group_weight)[group_index]
        background_only_rmse = float(np.sqrt(np.mean((background_only-observed)**2)))
        initial.extend([0.0] * (3 * len(group_names)))
        lower.extend([-np.inf] * (3 * len(group_names)))
        upper.extend([np.inf] * (3 * len(group_names)))
        names.extend(f"background_{group}.{axis}_nt" for group in group_names for axis in ("bx", "by", "bz"))
    if observed.size <= len(initial):
        raise ValueError("there must be more observed components than fitted parameters")
    library = library or build_fortran("orbital")

    def unpack(parameters):
        sources = deepcopy(initial_sources)
        for value, (index, key) in zip(parameters[:source_parameter_count], indices):
            sources[index][key] = float(value)
        return sources

    @lru_cache(maxsize=2)
    def crustal_field(parameters):
        return predict_field(xyz, unpack(parameters), radius_km=radius, quadrature_order=quadrature_order, library=library)

    def evaluate(parameters):
        sources = unpack(parameters)
        crustal = crustal_field(tuple(parameters[:source_parameter_count]))
        background = np.zeros_like(crustal)
        if fit_background:
            background = parameters[source_parameter_count:].reshape(-1, 3)[group_index]
        return sources, crustal, background

    def residual(parameters):
        sources, crustal, background = evaluate(parameters)
        values = ((crustal + background - observed) / sigma).ravel()
        penalties = [(sources[index][key] - mean) / std for index, prior in enumerate(priors)
                     for key, (mean, std) in prior.items()]
        return np.concatenate((values, penalties))

    initial_sources_copy, initial_prediction, initial_background = evaluate(np.asarray(initial))
    initial_rmse = float(np.sqrt(np.mean((initial_prediction + initial_background - observed)**2)))
    solution = least_squares(
        residual, initial, bounds=(lower, upper), method="trf", jac="3-point", x_scale="jac",
        loss=loss, f_scale=f_scale, max_nfev=max_nfev,
    )
    sources, crustal, background = evaluate(solution.x)
    prediction = crustal + background
    residual_nt = prediction - observed
    singular_values = np.linalg.svd(solution.jac, compute_uv=False)
    tolerance = max(solution.jac.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    return {
        "sources": sources, "predicted_nt": prediction, "crustal_nt": crustal,
        "background_nt": background, "residual_nt": residual_nt,
        "background_by_group": dict(zip(group_names, solution.x[source_parameter_count:].reshape(-1, 3).tolist()))
            if fit_background else {},
        "success": bool(solution.success), "message": str(solution.message), "nfev": int(solution.nfev),
        "cost": float(solution.cost), "initial_rmse_nt": initial_rmse,
        "background_only_rmse_nt": background_only_rmse,
        "rmse_nt": float(np.sqrt(np.mean(residual_nt**2))),
        "component_rmse_nt": np.sqrt(np.mean(residual_nt**2, axis=0)).tolist(),
        "jacobian_rank": rank, "singular_values": singular_values.tolist(),
        "parameter_names": names, "active_bounds": solution.active_mask.tolist(),
        "radius_km": radius, "quadrature_order": quadrature_order,
        "loss": loss, "f_scale": f_scale, "bounds": bounds, "priors": priors,
    }


def save_fit(result, observations, directory=None):
    """write source parameters, observation-aligned predictions, and fit metadata"""
    directory = (artifact_path(observation_input_path(observations)).parent if directory is None
                 else Path(directory).expanduser().resolve())
    directory.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(result["sources"]).to_csv(directory / "sources.csv", index=False)
    if "table" in observations:
        table = observations["table"].copy()
    else:
        table = pd.DataFrame(observations["xyz_km"], columns=["x_km", "y_km", "z_km"])
    for index, key in enumerate(("bx", "by", "bz")):
        table[f"observed_{key}_nt"] = observations["field_nt"][:, index]
        for name in ("predicted", "crustal", "background", "residual"):
            table[f"{name}_{key}_nt"] = result[f"{name}_nt"][:, index]
    table["observed_btot_nt"] = np.linalg.norm(observations["field_nt"], axis=1)
    table["predicted_btot_nt"] = np.linalg.norm(result["predicted_nt"], axis=1)
    table.to_csv(directory / "predictions.csv", index=False)
    metadata = {key: value for key, value in result.items() if not isinstance(value, np.ndarray)}
    metadata["data_report"] = observations.get("report", [])
    metadata["frame"] = observations.get("frame", "planet_fixed")
    altitude = np.linalg.norm(observations["xyz_km"], axis=1) - result["radius_km"]
    metadata["altitude_range_km"] = [float(altitude.min()), float(altitude.max())]
    (directory / "fit.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return directory
