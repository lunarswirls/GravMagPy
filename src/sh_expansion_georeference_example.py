#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Imports:
import numpy as np
import pyshtools as pysh
import rasterio
from rasterio.transform import from_origin
from rasterio.crs import CRS
import src.utils.wkt_defs as wkt_defs

"""
@author lunarswirls
"""
# Parameters
input_file = '/Users/wallecr1/Downloads/GRGM1200L_BOUGUER.sha'  # Replace with your file
output_tiff = '/Users/wallecr1/Downloads/bouguer_filtered_deg6_600.tif'
min_degree = 6
max_degree = 600

print("Loading spherical harmonic coefficients...")
# Load SH coeffs
coeffs = pysh.SHGravCoeffs.from_file(input_file)

print(f"Filtering coefficients: keeping degrees {min_degree} to {max_degree}...")
# Filter degrees
coeffs_filtered = coeffs.pad(max_degree)
coeffs_filtered.coeffs[:, :min_degree, :] = 0

print("Expanding spherical harmonics to spatial grid...")
# Expand to grid
grid = coeffs_filtered.expand(lmax=max_degree)
data = grid.total.data

# Get grid dimensions and resolution
nlat, nlon = data.shape
lat_res = 180.0 / nlat
lon_res = 360.0 / nlon

# Flip latitude to have north at the top
data = np.flipud(data)

# Set raster transform
transform = from_origin(west=0, north=90, xsize=lon_res, ysize=lat_res)

print("Georeferencing data...")

moon_wkt = wkt_defs.WKT_dict["GCS_Moon_2000"]

moon_crs = CRS.from_wkt(wkt=moon_wkt)

# Same as GeoTIFF
with rasterio.open(
    output_tiff,
    'w',
    driver='GTiff',
    height=data.shape[0],
    width=data.shape[1],
    count=1,
    dtype='float32',
    crs=moon_crs,
    transform=transform,
) as dst:
    dst.write(data.astype('float32'), 1)

print(f"GeoTIFF saved as {output_tiff}")
