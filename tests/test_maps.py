"""georeferencing, map alignment, and equivalent field continuation checks"""

import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
import pandas as pd

from gravmagpy import equivalent_field_grid, predict_equivalent_field, spherical_to_cartesian


def dipoles():
    return pd.DataFrame({"lat_deg": [0.0], "lon_deg": [0.0], "radius_m": [1727400.0],
                         "mx_am2": [1e10], "my_am2": [-2e10], "mz_am2": [3e10], "layer_id": [1]})


class equivalent_mapping_tests(unittest.TestCase):
    def test_volume_backend_does_not_require_new_dipole_entrypoint(self):
        from gravmagpy.forward import load_backend
        backend = SimpleNamespace(orbital_field=lambda: None)
        with patch("gravmagpy.forward.ctypes.CDLL", return_value=backend):
            self.assertIs(load_backend("older-volume-library"), backend)
        self.assertFalse(hasattr(backend, "dipole_field"))
        load_backend.cache_clear()

    def test_fortran_mapping_matches_independent_dipole_and_altitude_decay(self):
        table = dipoles()
        xyz = np.array([[1767.4, 0, 0], [1797.4, 0, 0]])
        actual = predict_equivalent_field(xyz, table)
        moment = table[["mx_am2", "my_am2", "mz_am2"]].to_numpy()[0]
        distances = np.array([40000.0, 70000.0])
        expected = 100*np.array([2*moment[0], -moment[1], -moment[2]]) / distances[:, None]**3
        np.testing.assert_allclose(actual, expected, rtol=1e-12)
        self.assertGreater(np.linalg.norm(actual[0]), np.linalg.norm(actual[1]))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "dipoles.csv"
            table.to_csv(path, index=False)
            np.testing.assert_allclose(predict_equivalent_field(xyz, path), expected)

    def test_map_pixel_centers_north_up_and_multilayer_sum(self):
        table = dipoles()
        both = pd.concat([table, table.assign(layer_id=2)], ignore_index=True)
        grid = equivalent_field_grid(both, bounds=(-2, -1, 2, 1), shape=(2, 4), altitude_km=30)
        np.testing.assert_allclose(grid["latitude_deg"], [0.5, -0.5])
        np.testing.assert_allclose(grid["longitude_deg"], [-1.5, -0.5, 0.5, 1.5])
        single = equivalent_field_grid(table, bounds=(-2, -1, 2, 1), shape=(2, 4), altitude_km=30)
        np.testing.assert_allclose(grid["field"], 2*single["field"])
        np.testing.assert_allclose(grid["btot_nt"], np.linalg.norm(grid["field"], axis=-1))
        with self.assertRaises(ValueError):
            equivalent_field_grid(table, bounds=(-2, -1, 2, 1), altitude_km=-1)
        with self.assertRaises(RuntimeError):
            predict_equivalent_field([[1737.4, 0, 0]], table.assign(radius_m=1737400.0))


@unittest.skipUnless(importlib.util.find_spec("rasterio"), "rasterio is an optional maps dependency")
class geotiff_mapping_tests(unittest.TestCase):
    def write_raster(self, directory, values, *, transform=None, crs=None, name="map.tif", nodata=-9999):
        import rasterio
        from rasterio.transform import from_bounds
        from gravmagpy.maps import geographic_crs
        path = Path(directory) / name
        values = np.asarray(values, dtype="float32")
        with rasterio.open(path, "w", driver="GTiff", count=1, height=values.shape[0], width=values.shape[1],
                           dtype="float32", transform=transform or from_bounds(-2, -2, 2, 2, values.shape[1], values.shape[0]),
                           crs=geographic_crs() if crs is None else crs, nodata=nodata) as target:
            target.write(values, 1)
        return path

    def test_pixel_alignment_nodata_and_scale_metadata(self):
        import rasterio
        from gravmagpy.maps import read_geotiff
        values = np.arange(16).reshape(4, 4).astype(float)
        values[0, 0] = -9999
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, values)
            with rasterio.open(path, "r+") as target:
                target.scales = (0.1,)
                target.offsets = (2.0,)
                mask = np.full((4, 4), 255, dtype="uint8")
                mask[1, 1] = 0
                target.write_mask(mask)
            result = read_geotiff(path, bounds=(-2, -2, 2, 2), shape=(4, 4), resampling="nearest")
            self.assertTrue(result["data"].mask[0, 0])
            self.assertTrue(result["data"].mask[1, 1])
            np.testing.assert_allclose(result["data"][2:], values[2:]*0.1+2)
            np.testing.assert_allclose(result["latitude_deg"], [1.5, 0.5, -0.5, -1.5])
            self.assertEqual(result["transform"].e, -1)

    def test_projected_lunar_raster_matches_geographic_grid(self):
        from rasterio.crs import CRS
        from rasterio.transform import from_bounds
        from gravmagpy.maps import read_geotiff
        radius = 1737400.0
        edge = np.deg2rad(2)*radius
        values = np.arange(16).reshape(4, 4)
        crs = CRS.from_string(f"+proj=eqc +R={radius} +lat_ts=0 +lon_0=0 +units=m +no_defs")
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, values, crs=crs, transform=from_bounds(-edge, -edge, edge, edge, 4, 4))
            result = read_geotiff(path, bounds=(-2, -2, 2, 2), shape=(4, 4), resampling="nearest")
            np.testing.assert_allclose(result["data"], values)

    def test_zero_to_360_and_antimeridian(self):
        from rasterio.transform import from_bounds
        from gravmagpy.maps import read_geotiff
        values = np.tile(np.arange(8), (4, 1))
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, values, transform=from_bounds(0, -20, 360, 20, 8, 4))
            result = read_geotiff(path, bounds=(-100, -10, -80, 10), shape=(4, 4), resampling="nearest")
            np.testing.assert_allclose(result["data"], np.tile([5, 5, 6, 6], (4, 1)))
            result = read_geotiff(path, bounds=(-10, -10, 10, 10), shape=(4, 4), resampling="nearest")
            self.assertFalse(np.ma.getmaskarray(result["data"]).any())
            np.testing.assert_allclose(result["data"], np.tile([7, 7, 0, 0], (4, 1)))
            path = self.write_raster(directory, values, transform=from_bounds(-180, -20, 180, 20, 8, 4), name="seam.tif")
            result = read_geotiff(path, bounds=(170, -10, -170, 10), shape=(4, 4), resampling="nearest")
            np.testing.assert_allclose(result["data"], np.tile([7, 7, 0, 0], (4, 1)))
            np.testing.assert_allclose(result["longitude_deg"], [172.5, 177.5, 182.5, 187.5])

    def test_rotated_raster_uses_affine_transform(self):
        from rasterio.transform import Affine, rowcol
        from gravmagpy.maps import read_geotiff
        transform = Affine.translation(-2, 2)*Affine.rotation(15)*Affine.scale(1, -1)
        values = np.arange(16).reshape(4, 4)
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, values, transform=transform)
            result = read_geotiff(path, bounds=(-3, -3, 3, 3), shape=(12, 12), resampling="nearest")
            lon, lat = np.meshgrid(result["longitude_deg"], result["latitude_deg"])
            rows, cols = np.array(rowcol(transform, lon.ravel(), lat.ravel()))
            inside = (rows >= 0) & (rows < 4) & (cols >= 0) & (cols < 4)
            expected = values[rows[inside], cols[inside]]
            np.testing.assert_allclose(result["data"].ravel()[inside], expected)

    def test_missing_wrong_crs_no_overlap_and_invalid_settings(self):
        import rasterio
        from rasterio.transform import from_bounds
        from gravmagpy.maps import read_geotiff, geographic_crs
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, np.ones((4, 4)), crs="EPSG:4326")
            with self.assertRaises(ValueError):
                read_geotiff(path, bounds=(-2, -2, 2, 2))
            result = read_geotiff(path, bounds=(-2, -2, 2, 2), shape=(4, 4), source_crs=geographic_crs())
            np.testing.assert_allclose(result["data"], 1)
            path = Path(directory) / "missing.tif"
            with rasterio.open(path, "w", driver="GTiff", count=1, height=4, width=4, dtype="float32",
                               transform=from_bounds(-2, -2, 2, 2, 4, 4)) as target:
                target.write(np.ones((4, 4), dtype="float32"), 1)
            with self.assertRaises(ValueError):
                read_geotiff(path, bounds=(-2, -2, 2, 2))
            with self.assertRaises(ValueError):
                read_geotiff(path, bounds=(30, 30, 40, 40), source_crs=geographic_crs())
            with self.assertRaises(ValueError):
                read_geotiff(path, bounds=(-2, -2, 2, 2), band=2, source_crs=geographic_crs())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib"), "matplotlib is an optional maps dependency")
    def test_plot_panels_with_contours_and_saved_sources(self):
        from gravmagpy.maps import plot_equivalent_maps
        with tempfile.TemporaryDirectory() as directory:
            path = self.write_raster(directory, np.arange(16).reshape(4, 4))
            source_path = Path(directory) / "dipoles.csv"
            dipoles().to_csv(source_path, index=False)
            output = Path(directory) / "comparison.png"
            result = plot_equivalent_maps(source_path, [{"path": path, "title": "synthetic spectral map", "unit": "ratio"}],
                                          output, bounds=(-2, -2, 2, 2), shape=(15, 15), components=("bx", "btot"), dpi=60)
            self.assertEqual(result, output.resolve())
            self.assertGreater(output.stat().st_size, 1000)
            with self.assertRaises(ValueError):
                plot_equivalent_maps(source_path, [{"path": path}], path, bounds=(-2, -2, 2, 2))
