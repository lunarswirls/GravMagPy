#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Imports:
import numpy as np
import rasterio
from scipy.optimize import minimize
from numba import njit, prange

"""
@author lunarswirls
"""
# Define constants
MU_0 = 4 * np.pi * 1e-7  # Permeability of free space (T·m/A)


def load_geotiff(filename):
    """
    Load a GeoTIFF file and extract data and georeferencing information.
    
    Parameters:
    filename: str
        Path to the GeoTIFF file.
    
    Returns:
    data: np.array
        2D array of magnetic field intensity.
    transform: Affine
        Geotransformation matrix.
    crs: dict
        Coordinate reference system.
    """
    with rasterio.open(filename) as src:
        data = src.read(1)  # Read first band (magnetic field)
        transform = src.transform
        crs = src.crs  # Coordinate reference system
        return data, transform, crs


def latlon_to_indices(lat, lon, transform):
    """
    Convert latitude and longitude to row and column indices in the GeoTIFF.
    
    Parameters:
    lat, lon: float
        Latitude and Longitude of the desired point.
    transform: Affine
        Geotransformation matrix.
    
    Returns:
    row, col: int
        Pixel row and column indices.
    """
    col, row = ~transform * (lon, lat)  # Convert lat/lon to pixel indices
    return int(row), int(col)


def extract_subset(data_array, tf, lat_range, lon_range):
    """
    Extract a subset of magnetic field data within a given lat/lon range.
    
    Parameters:
    data_array: np.array
        Full magnetic field dataset.
    tf: Affine
        Geotransformation matrix.
    lat_range, lon_range: tuple
        Min and max values of latitude and longitude.
    
    Returns:
    subset: np.array
        Extracted subset of magnetic field data.
    coords: list of tuples
        List of (latitude, longitude) corresponding to the extracted data.
    """
    row_min, col_min = latlon_to_indices(lat_range[0], lon_range[0], tf)
    row_max, col_max = latlon_to_indices(lat_range[1], lon_range[1], tf)

    subset = data_array[row_min:row_max, col_min:col_max]
    
    # Compute corresponding lat/lon coordinates
    coords = []
    for row in range(row_min, row_max):
        for col in range(col_min, col_max):
            lon, lat = tf * (col, row)
            coords.append((lat, lon))
    
    return subset, np.array(coords)


@njit(fastmath=True, parallel=True)
def prism_magnetic_field_analytical(params, obs_points, n_prisms):
    """
    Compute the total magnetic field at observation points due to multiple magnetized prisms.
    Optimized using Numba for performance.
    """
    B_field = np.zeros_like(obs_points)

    for n in prange(n_prisms):
        offset = n * 9
        x_min, x_max, y_min, y_max, z_min, z_max, Mx, My, Mz = params[offset:offset+9]

        for i in prange(len(obs_points)):
            x, y, z = obs_points[i]
            Bx, By, Bz = 0.0, 0.0, 0.0

            for dx in (x_min, x_max):
                for dy in (y_min, y_max):
                    for dz in (z_min, z_max):
                        rx, ry, rz = x - dx, y - dy, z - dz
                        R = np.sqrt(rx**2 + ry**2 + rz**2) + 1e-6

                        sign = (-1) ** ((dx == x_max) + (dy == y_max) + (dz == z_max))

                        Bx += sign * np.arctan2(ry * rz, rx * R)
                        By += sign * np.arctan2(rx * rz, ry * R)
                        Bz += sign * np.arctan2(rx * ry, rz * R)

            B_field[i, 0] += MU_0 * Mx * Bx
            B_field[i, 1] += MU_0 * My * By
            B_field[i, 2] += MU_0 * Mz * Bz

    return B_field


def fit_prism_sources(obs_points, B_obs, n_prisms):
    """
    Perform optimization to find the best-fit prism parameters.
    """
    initial_guess = []
    for _ in range(n_prisms):
        x_min, x_max = np.sort(np.random.uniform(-30, 30, 2))
        y_min, y_max = np.sort(np.random.uniform(-30, 30, 2))
        z_min, z_max = np.sort(np.random.uniform(-40, -5, 2))
        Mx, My, Mz = np.random.uniform(-2, 2, 3)
        initial_guess.extend([x_min, x_max, y_min, y_max, z_min, z_max, Mx, My, Mz])

    initial_guess = np.array(initial_guess)

    bounds = [(-50, 50), (-50, 50), (-50, 50), (-50, 50), (-50, 0), (-50, 0), (-5, 5), (-5, 5), (-5, 5)] * n_prisms
    
    result = minimize(lambda params: np.sum((B_obs - prism_magnetic_field_analytical(params, obs_points, num_prisms))**2),
                      initial_guess, bounds=bounds, method='L-BFGS-B')
    return result.x  # Optimized parameters


# Main execution
fname = "/Users/danywaller/Data/LP_modeled_surface_total_field_2ppd.tiff"  # Example GeoTIFF file
lt_range = (-10, 10)  # Example latitude range
ln_range = (-20, 20)  # Example longitude range

# Load GeoTIFF
data, transform, crs = load_geotiff(fname)

# Extract a subset of magnetic field data
subset, coords = extract_subset(data, transform, lt_range, ln_range)

# Convert to Cartesian coordinates for modeling
observation_points = np.array([(lat, lon, 30) for lat, lon in coords])  # Assume 30 km altitude

# Generate synthetic observed magnetic field data
num_prisms = 5
B_observed = np.random.normal(0, 0.1, (len(observation_points), 3))  # Placeholder data

# Fit prism sources
best_fit_params = fit_prism_sources(observation_points, B_observed, num_prisms)

print("\nBest-Fit Parameters:\n", best_fit_params.reshape(num_prisms, 9))
