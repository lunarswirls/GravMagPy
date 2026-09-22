"""published source parameters, moments, frames, and numerical approximations"""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from gravmagpy import read_model_input, run_sphere_model, spherical_to_cartesian, write_model_input
from examples.reiner_gamma_hemingway_garrick_bethell_2012_test import model as hemingway
from examples.reiner_gamma_chaffee_2025_test import model as chaffee


def point_dipole_field(latitude_deg, longitude_deg, altitude_km):
    """independent analytic point-dipole reference in the published local frames"""
    lat, lon = np.meshgrid(latitude_deg, longitude_deg, indexing="ij")
    xyz = spherical_to_cartesian(lat, lon, 1737.4+altitude_km)*1000
    field = np.zeros_like(xyz)
    for row in hemingway.load_dipoles():
        position = spherical_to_cartesian(row["latitude_deg"], row["longitude_deg"], 1737.4-row["depth_km"])*1000
        moment = row["moment_a_m2"]*hemingway.dipole_direction(row)
        delta = xyz-position
        distance = np.linalg.norm(delta, axis=-1)
        field += 100*(3*np.sum(delta*moment, axis=-1)[..., None]*delta/distance[..., None]**5
                      - moment/distance[..., None]**3)
    return field


def ellipsoid_volume_field(latitude_deg, longitude_deg, altitude_km, order=20):
    """independently integrate the mapped ellipsoid volumes using dipole kernels"""
    parameters = chaffee.load_parameters()
    lat, lon = np.meshgrid(latitude_deg, longitude_deg, indexing="ij")
    observers = spherical_to_cartesian(lat.ravel(), lon.ravel(), 1737.4+altitude_km)*1000
    field = np.zeros_like(observers)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    z, radial, theta = np.meshgrid(nodes, (nodes+1)/2, np.linspace(0, 2*np.pi, 128, endpoint=False), indexing="ij")
    wz, wr, wt = np.meshgrid(weights, weights/2, np.full(128, 2*np.pi/128), indexing="ij")
    for body in parameters["bodies"]:
        a, b, c = body["semiaxes_km"]
        major = a*np.sqrt(1-z*z)*radial*np.cos(theta)
        minor = b*np.sqrt(1-z*z)*radial*np.sin(theta)
        strike = np.deg2rad(body["strike_deg"])
        north = body["north_km"]+major*np.cos(strike)-minor*np.sin(strike)
        east = body["east_km"]+major*np.sin(strike)+minor*np.cos(strike)
        coordinates = chaffee.surface_coordinates(north.ravel(), east.ravel(), origin_latitude_deg=7.5,
                                                  origin_longitude_deg=-59, radius_km=1737.4)
        positions = spherical_to_cartesian(coordinates[:, 0], coordinates[:, 1],
                                           1737.4-body["depth_center_km"]-c*z.ravel())*1000
        volumes = (a*b*c*(1-z*z)*radial*wz*wr*wt).ravel()*1e9
        moments = volumes[:, None]*chaffee.magnetization_vector(parameters)
        for index, observer in enumerate(observers):
            delta = observer-positions
            distance = np.linalg.norm(delta, axis=1)
            field[index] += 100*np.sum(3*np.sum(delta*moments, axis=1)[:, None]*delta/distance[:, None]**5
                                       - moments/distance[:, None]**3, axis=0)
    return field.reshape(*lat.shape, 3)


class reiner_gamma_tests(unittest.TestCase):
    def test_hemingway_published_moments_and_local_directions(self):
        rows = hemingway.load_dipoles()
        sources = hemingway.build_model()["sources"]
        self.assertEqual(len(rows), 55)
        self.assertAlmostEqual(sum(row["moment_a_m2"] for row in rows), 1.001e13)
        self.assertEqual((rows[-1]["latitude_deg"], rows[-1]["longitude_deg"]), (7.557, 301.749))
        for row, source in zip(rows, sources):
            lat_low, lat_high = np.deg2rad(source["latitude_deg"])
            lon_low, lon_high = np.deg2rad(source["longitude_deg"])
            top, bottom = (1737.4-np.asarray(source["depth_km"]))*1000
            volume = (lon_high-lon_low)*(np.sin(lat_high)-np.sin(lat_low))*(top**3-bottom**3)/3
            moment = volume*np.asarray(source["magnetization_a_m"])
            np.testing.assert_allclose(np.linalg.norm(moment), row["moment_a_m2"], rtol=2e-11)
            lat, lon = np.deg2rad([row["latitude_deg"], row["longitude_deg"]])
            basis = np.array([[-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
                              [-np.sin(lon), np.cos(lon), 0],
                              [-np.cos(lat)*np.cos(lon), -np.cos(lat)*np.sin(lon), -np.sin(lat)]])
            north, east, down = basis @ moment
            self.assertAlmostEqual(np.rad2deg(np.arctan2(down, np.hypot(north, east))), 2.0)
            self.assertAlmostEqual(np.rad2deg(np.arctan2(east, north)), -8.0)
            np.testing.assert_allclose(source["depth_km"], [4.8, 5.2])

    def test_hemingway_field_converges_to_published_point_dipoles(self):
        latitude = np.linspace(7.0, 7.8, 5)
        longitude = np.linspace(-59.7, -57.9, 7)
        for altitude in (0.0, 18.0):
            reference = point_dipole_field(latitude, longitude, altitude)
            model = hemingway.build_model(altitude_km=altitude, latitude_deg=latitude, longitude_deg=longitude)
            result = run_sphere_model(model, solver_options={"refine_factor": 2})
            relative_error = np.linalg.norm(result["field"]-reference)/np.linalg.norm(reference)
            self.assertLess(relative_error, 0.005)
        large = hemingway.build_model(altitude_km=0, latitude_deg=latitude, longitude_deg=longitude,
                                     box_width_deg=0.08, thickness_km=2)
        large_field = run_sphere_model(large, solver_options={"refine_factor": 2})["field"]
        small = hemingway.build_model(altitude_km=0, latitude_deg=latitude, longitude_deg=longitude)
        small_field = run_sphere_model(small, solver_options={"refine_factor": 2})["field"]
        reference = point_dipole_field(latitude, longitude, 0)
        self.assertLess(np.linalg.norm(small_field-reference), np.linalg.norm(large_field-reference))

    def test_chaffee_parameters_match_released_thermal_snapshot(self):
        parameters = chaffee.load_parameters()
        data = np.loadtxt(Path(chaffee.__file__).with_name("reiner-smalltube-buried-02.dat"), skiprows=7)
        hot = data[data[:, 2] >= parameters["thermal_threshold_c"]]
        width, height = np.ptp(hot[:, :2], axis=0)*parameters["thermal_cell_m"]/1000
        depth = -(hot[:, 1].min()*parameters["thermal_cell_m"]/1000+height/2)
        np.testing.assert_allclose([width, height, depth], [5, 1.5, 1.75])
        self.assertEqual([body["strike_deg"] for body in parameters["bodies"]], [70, 70, 40, 80])
        for body in parameters["bodies"]:
            np.testing.assert_allclose(body["semiaxes_km"][1:], [width, height])
            self.assertEqual(body["depth_center_km"], depth)

    def test_chaffee_slice_volumes_and_common_magnetization(self):
        parameters = chaffee.load_parameters()
        for layers in (8, 16):
            for body in parameters["bodies"]:
                sources = chaffee.ellipsoid_layers(body, parameters, layers=layers)
                volume = sum(chaffee.polygon_volume(source["vertices_lat_lon"], source["depth_km"]) for source in sources)
                expected = 4*np.pi*np.prod(body["semiaxes_km"])*1e9/3
                np.testing.assert_allclose(volume, expected, rtol=1e-9)
                self.assertAlmostEqual(sources[0]["depth_km"][0], 0.25)
                self.assertAlmostEqual(sources[-1]["depth_km"][1], 3.25)
                for source in sources:
                    np.testing.assert_allclose(np.linalg.norm(source["magnetization_a_m"]), 0.5, atol=1e-14)
                    np.testing.assert_allclose(source["magnetization_a_m"], chaffee.magnetization_vector(parameters))

    def test_polygon_volume_and_geographic_mapping(self):
        outline = [[6, -61], [6, -60], [7, -60], [7, -61]]
        expected = np.deg2rad(1)*(np.sin(np.deg2rad(7))-np.sin(np.deg2rad(6)))
        expected *= ((1737.4-1)**3-(1737.4-2)**3)*1e9/3
        np.testing.assert_allclose(chaffee.polygon_volume(outline, [1, 2]), expected, rtol=1e-11)
        coordinates = chaffee.surface_coordinates([0, 1737.4*np.pi/180], [0, 0],
                                                  origin_latitude_deg=7.5, origin_longitude_deg=-59, radius_km=1737.4)
        np.testing.assert_allclose(coordinates, [[7.5, -59], [8.5, -59]], atol=1e-12)

    def test_both_models_round_trip_through_input_cards(self):
        with tempfile.TemporaryDirectory() as directory:
            for index, model in enumerate((hemingway.build_model(), chaffee.build_model(layers=4, vertices=16))):
                path = write_model_input(model, Path(directory)/f"model_{index}.in")
                restored = read_model_input(path)
                self.assertEqual(len(restored["sources"]), len(model["sources"]))
                self.assertEqual(restored["grid"], model["grid"])
                for original, actual in zip(model["sources"], restored["sources"]):
                    np.testing.assert_allclose(actual["magnetization_a_m"], original["magnetization_a_m"], rtol=1e-9)
                    np.testing.assert_allclose(actual["depth_km"], original["depth_km"])

    def test_chaffee_field_resolution_at_orbital_height(self):
        grid = {"latitude_deg": np.linspace(6, 9, 4), "longitude_deg": np.linspace(-61, -57, 5)}
        model = chaffee.build_model(**grid, layers=16, vertices=64, mesh_spacing_km=0.75)
        finer = chaffee.build_model(**grid, layers=32, vertices=128, mesh_spacing_km=0.375)
        field = run_sphere_model(model, solver_options={"refine_factor": 1})["field"]
        reference = run_sphere_model(finer, solver_options={"refine_factor": 1})["field"]
        self.assertGreater(np.linalg.norm(field), 0)
        self.assertLess(np.linalg.norm(field-reference)/np.linalg.norm(reference), 0.03)
        volume_reference = ellipsoid_volume_field(**grid, altitude_km=30)
        self.assertLess(np.linalg.norm(field-volume_reference)/np.linalg.norm(volume_reference), 0.03)

    def test_chaffee_surface_field_resolution(self):
        grid = {"latitude_deg": np.linspace(7.2, 7.8, 4), "longitude_deg": np.linspace(-59.3, -58.7, 5),
                "altitude_km": 0}
        model = chaffee.build_model(**grid)
        finer = chaffee.build_model(**grid, layers=128, vertices=256, mesh_spacing_km=0.09375)
        field = run_sphere_model(model, solver_options={"refine_factor": 1})["field"]
        reference = run_sphere_model(finer, solver_options={"refine_factor": 1})["field"]
        self.assertLess(np.linalg.norm(field-reference)/np.linalg.norm(reference), 0.03)
