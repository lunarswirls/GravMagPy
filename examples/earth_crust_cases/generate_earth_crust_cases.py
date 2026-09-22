from pathlib import Path
import math


output_dir = Path(__file__).resolve().parent
earth_radius_km = 6371.2


def wrap180(lon_deg):
    return ((lon_deg + 180.0) % 360.0) - 180.0


def local_geophysical_to_global_angles(lat_deg, lon_deg, inc_local_deg, dec_local_deg):
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    inc = math.radians(inc_local_deg)
    dec = math.radians(dec_local_deg)

    north = math.cos(inc) * math.cos(dec)
    east = math.cos(inc) * math.sin(dec)
    down = math.sin(inc)

    east_axis = (-math.sin(lon), math.cos(lon), 0.0)
    north_axis = (
        -math.sin(lat) * math.cos(lon),
        -math.sin(lat) * math.sin(lon),
        math.cos(lat),
    )
    up_axis = (
        math.cos(lat) * math.cos(lon),
        math.cos(lat) * math.sin(lon),
        math.sin(lat),
    )

    mx = north * north_axis[0] + east * east_axis[0] - down * up_axis[0]
    my = north * north_axis[1] + east * east_axis[1] - down * up_axis[1]
    mz = north * north_axis[2] + east * east_axis[2] - down * up_axis[2]

    inc_global_deg = math.degrees(math.atan2(mz, math.hypot(mx, my)))
    dec_global_deg = math.degrees(math.atan2(my, mx))
    return inc_global_deg, dec_global_deg


def axial_dipole_local_inclination(lat_deg):
    lat = math.radians(lat_deg)
    return math.degrees(math.atan2(2.0 * math.sin(lat), math.cos(lat)))


def great_circle_distance_km(lat1_deg, lon1_deg, lat2_deg, lon2_deg):
    lat1 = math.radians(lat1_deg)
    lat2 = math.radians(lat2_deg)
    dlat = lat2 - lat1
    dlon = math.radians(lon2_deg - lon1_deg)
    hav = math.sin(0.5 * dlat) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(0.5 * dlon) ** 2
    return 2.0 * earth_radius_km * math.asin(min(1.0, math.sqrt(hav)))


def write_case(path, bodies, grid):
    with path.open("w") as handle:
        for index, body in enumerate(bodies, start=1):
            handle.write(f"{body['title']} source {index}\n")
            handle.write(
                f"{grid['lat0']:.6f} {grid['lon0']:.6f} "
                f"{grid['dlat']:.6f} {grid['dlon']:.6f} {grid['elevation_km']:.6f} "
                f"{grid['nlat']} {grid['nlon']}\n"
            )
            handle.write("1 1 1 1\n")
            handle.write("0.1 1.50 -2.0 2.0 0.080\n")
            handle.write(
                f"2 1 {body['magnetization_apm']:.9e} "
                f"{body['inc_global_deg']:.6f} {body['dec_global_deg']:.6f}\n"
            )
            handle.write("0 1 1 0.0 0.0\n")
            handle.write(
                f"{body['lat_max']:.6f} {body['lat_min']:.6f} "
                f"{body['lon_max']:.6f} {body['lon_min']:.6f} "
                f"{body['depth_top_km']:.6f} {body['depth_bottom_km']:.6f}\n\n"
            )


def build_bangui_case():
    config = {
        "center_lat_deg": 4.37,
        "center_lon_deg": 18.56,
        "radius_km": 400.0,
        "spacing_deg": 0.5,
        "depth_top_km": 3.0,
        "depth_bottom_km": 7.5,
        "magnetization_apm": 10.0,
        "inc_local_deg": 25.0,
        "dec_local_deg": -18.0,
    }
    grid = {
        "lat0": -5.0,
        "lon0": 8.0,
        "dlat": 0.25,
        "dlon": 0.25,
        "elevation_km": 4.0,
        "nlat": 81,
        "nlon": 89,
    }
    bodies = []
    half = 0.5 * config["spacing_deg"]
    for ilat in range(41):
        lat = -5.0 + ilat * config["spacing_deg"]
        for ilon in range(45):
            lon = 8.0 + ilon * config["spacing_deg"]
            distance_km = great_circle_distance_km(
                lat,
                lon,
                config["center_lat_deg"],
                config["center_lon_deg"],
            )
            if distance_km > config["radius_km"]:
                continue
            inc_global_deg, dec_global_deg = local_geophysical_to_global_angles(
                lat,
                lon,
                config["inc_local_deg"],
                config["dec_local_deg"],
            )
            bodies.append(
                {
                    "title": "bangui published disc proxy",
                    "lat_max": lat + half,
                    "lat_min": lat - half,
                    "lon_max": lon + half,
                    "lon_min": lon - half,
                    "depth_top_km": config["depth_top_km"],
                    "depth_bottom_km": config["depth_bottom_km"],
                    "magnetization_apm": config["magnetization_apm"],
                    "inc_global_deg": inc_global_deg,
                    "dec_global_deg": dec_global_deg,
                }
            )
    return bodies, grid


def global_class_index(lat_deg, lon_deg, layer_id):
    lat_band = int((lat_deg + 90.0) // 18.0)
    lon_band = int((lon_deg + 180.0) // 36.0)
    return (lat_band + 3 * lon_band + 2 * layer_id) % 10


def build_meyer_global_case():
    grid = {
        "lat0": -90.0,
        "lon0": -180.0,
        "dlat": 5.0,
        "dlon": 5.0,
        "elevation_km": 400.0,
        "nlat": 37,
        "nlon": 73,
    }
    layer_limits_km = [(0.0, 20.0), (20.0, 40.0)]
    class_magnetization_apm = [0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.48, 0.65, 0.85, 1.10]
    bodies = []
    for layer_id, limits in enumerate(layer_limits_km):
        for ilat in range(90):
            lat = -89.0 + 2.0 * ilat
            inc_local_deg = axial_dipole_local_inclination(lat)
            for ilon in range(180):
                lon = -179.0 + 2.0 * ilon
                class_index = global_class_index(lat, lon, layer_id)
                inc_global_deg, dec_global_deg = local_geophysical_to_global_angles(
                    lat,
                    lon,
                    inc_local_deg,
                    0.0,
                )
                bodies.append(
                    {
                        "title": f"meyer geometry proxy layer {layer_id + 1} class {class_index + 1}",
                        "lat_max": min(90.0, lat + 1.0),
                        "lat_min": max(-90.0, lat - 1.0),
                        "lon_max": lon + 1.0,
                        "lon_min": lon - 1.0,
                        "depth_top_km": limits[0],
                        "depth_bottom_km": limits[1],
                        "magnetization_apm": class_magnetization_apm[class_index],
                        "inc_global_deg": inc_global_deg,
                        "dec_global_deg": dec_global_deg,
                    }
                )
    return bodies, grid


def north_america_magnetization_apm(lat_deg, lon_deg):
    if 45.0 <= lat_deg <= 62.0 and -108.0 <= lon_deg <= -76.0:
        return 2.50
    if lat_deg >= 48.0 and -125.0 <= lon_deg <= -55.0:
        return 1.60
    if lat_deg <= 50.0 and lon_deg <= -105.0:
        return 0.25
    if lat_deg <= 50.0 and lon_deg >= -85.0:
        return 1.00
    return 0.75


def build_north_america_case():
    grid = {
        "lat0": 10.0,
        "lon0": -175.0,
        "dlat": 2.0,
        "dlon": 2.0,
        "elevation_km": 400.0,
        "nlat": 36,
        "nlon": 66,
    }
    bodies = []
    lat = 16.0
    while lat <= 74.0 + 1.0e-9:
        lon_spacing = 2.0 / max(0.26, math.cos(math.radians(lat)))
        nlon = max(1, int(round(120.0 / lon_spacing)))
        lon_spacing = 120.0 / nlon
        for ilon in range(nlon):
            lon = -170.0 + (ilon + 0.5) * lon_spacing
            inc_local_deg = axial_dipole_local_inclination(lat)
            inc_global_deg, dec_global_deg = local_geophysical_to_global_angles(
                lat,
                lon,
                inc_local_deg,
                0.0,
            )
            bodies.append(
                {
                    "title": "mayhew north america equal area proxy",
                    "lat_max": lat + 1.0,
                    "lat_min": lat - 1.0,
                    "lon_max": lon + 0.5 * lon_spacing,
                    "lon_min": lon - 0.5 * lon_spacing,
                    "depth_top_km": 0.0,
                    "depth_bottom_km": 40.0,
                    "magnetization_apm": north_america_magnetization_apm(lat, lon),
                    "inc_global_deg": inc_global_deg,
                    "dec_global_deg": dec_global_deg,
                }
            )
        lat += 2.0
    return bodies, grid


cases = {
    "bangui_equivalent_dipoles.in": build_bangui_case(),
    "meyer_global_32400_dipoles.in": build_meyer_global_case(),
    "mayhew_north_america_equivalent_layer.in": build_north_america_case(),
}

for filename, case in cases.items():
    case_bodies, case_grid = case
    write_case(output_dir / filename, case_bodies, case_grid)
    print(f"wrote {filename}: {len(case_bodies)} sources")
