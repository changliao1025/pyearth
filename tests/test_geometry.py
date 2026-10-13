"""Unit tests for GIS geometry algorithms, polygon areas, angles, and simplification."""

import numpy as np
import pytest
from osgeo import ogr

from pyearth.gis.geometry import (
    calculate_angle_between_point,
    calculate_angle_between_point_normal,
    calculate_angle_between_vectors_degrees,
    calculate_distance_based_on_longitude_latitude,
    calculate_distance_to_line,
    calculate_distance_to_plane,
    calculate_intersect_on_great_circle,
    calculate_polygon_area,
    calculate_spherical_triangle_area,
    check_counter_clockwise,
    clean_geometry,
    convert_180_to_360,
    convert_360_to_180,
    create_box_from_longitude_latitude,
    douglas_peucker_geodetic,
    extract_unique_vertices_and_connectivity,
    find_great_circle_intersection_with_meridian,
    find_nearest_resolution,
    get_output_geometry_type,
    get_polygon_orientation,
    haversine,
    is_convex_polygon,
    project_point_onto_plane,
    visvalingam_whyatt_geodetic,
)


class TestAngles:
    def test_angle_between_vectors(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([0.0, 1.0, 0.0])
        assert pytest.approx(calculate_angle_between_vectors_degrees(v1, v2)) == 90.0

        # Collinear / opposite
        v3 = np.array([-1.0, 0.0, 0.0])
        assert pytest.approx(calculate_angle_between_vectors_degrees(v1, v3)) == 180.0

        # Zero vector
        assert calculate_angle_between_vectors_degrees(v1, np.zeros(3)) == 0.0

    def test_angle_between_points_on_sphere(self):
        # Angle at North Pole between chord vectors to (0,0) and (90,0) in 3D Euclidean space
        ang = calculate_angle_between_point(
            0.0, 0.0, 0.0, 90.0, 90.0, 0.0
        )
        assert pytest.approx(ang) == 60.0

    def test_angle_between_point_normal(self):
        ang = calculate_angle_between_point_normal(
            0.0, 0.0, 0.0, 90.0, 90.0, 0.0
        )
        assert 0.0 <= ang < 360.0


class TestDistances:
    def test_haversine_distance(self):
        # Quarter circumference on equator (0,0) to (90,0)
        d = calculate_distance_based_on_longitude_latitude(0.0, 0.0, 90.0, 0.0)
        assert pytest.approx(d, rel=1e-3) == 10018754.0

    def test_distance_to_line(self):
        # Line along equator from (0,0) to (10,0)
        # Query point at (5, 1) -> distance approx 1 degree in meters ~ 111 km
        dist = calculate_distance_to_line(0.0, 0.0, 5.0, 1.0, 10.0, 0.0)
        assert 100000.0 < dist < 120000.0

    def test_distance_to_great_circle_plane(self):
        # Plane defined by equator: (0,0) and (90,0) and Earth center
        # Point on equator: distance should be 0.0
        dist_on_plane = calculate_distance_to_plane(0.0, 0.0, 45.0, 0.0, 90.0, 0.0)
        assert pytest.approx(dist_on_plane, abs=1e-10) == 0.0

        # Point at North Pole (0, 90): distance to equatorial plane should be 1.0 (unit sphere)
        dist_pole = calculate_distance_to_plane(0.0, 0.0, 0.0, 90.0, 90.0, 0.0)
        assert pytest.approx(dist_pole, abs=1e-10) == 1.0


class TestGreatCircleOperations:
    def test_project_point_onto_plane(self):
        p = np.array([1.0, 2.0, 3.0])
        normal = np.array([0.0, 0.0, 1.0])
        proj = project_point_onto_plane(p, normal)
        np.testing.assert_allclose(proj, [1.0, 2.0, 0.0])

    def test_calculate_intersect_on_great_circle(self):
        # Closest point on equator between (0,0) and (90,0) to North Pole (0,90)
        lon, lat = calculate_intersect_on_great_circle(0.0, 0.0, 0.0, 90.0, 90.0, 0.0)
        assert pytest.approx(lat, abs=1e-5) == 0.0

    def test_meridian_intersection(self):
        # Great circle along equator crosses 45E at 0N
        lon, lat = find_great_circle_intersection_with_meridian(0.0, 0.0, 90.0, 0.0, 45.0)
        assert pytest.approx(lon) == 45.0
        assert pytest.approx(lat, abs=1e-5) == 0.0


class TestPolygonAreas:
    def test_haversine_func(self):
        assert haversine(0.0) == 0.0
        assert pytest.approx(haversine(np.pi)) == 1.0
        assert pytest.approx(haversine(np.pi / 2.0)) == 0.5

    def test_spherical_triangle_area(self):
        # Octant of sphere in radians: lon=[0, pi/2, 0], lat=[0, 0, pi/2]
        lons_rad = [0.0, np.pi / 2.0, 0.0]
        lats_rad = [0.0, 0.0, np.pi / 2.0]
        area = calculate_spherical_triangle_area(lons_rad, lats_rad, iFlag_radian=True)
        assert pytest.approx(area, rel=1e-4) == np.pi / 2.0

    def test_polygon_area_algorithms(self):
        # 1-degree square at equator
        lon = [0.0, 1.0, 1.0, 0.0, 0.0]
        lat = [0.0, 0.0, 1.0, 1.0, 0.0]

        area_karney = calculate_polygon_area(lon, lat, iFlag_algorithm=2)
        area_green = calculate_polygon_area(lon, lat, iFlag_algorithm=0)
        area_lhuilier = calculate_polygon_area(lon, lat, iFlag_algorithm=1)

        assert area_karney > 0.0
        assert pytest.approx(area_green, rel=0.1) == area_karney
        assert pytest.approx(area_lhuilier, rel=1e-5) == area_karney


class TestPolygonChecksAndClean:
    def test_convex_polygon(self):
        square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
        assert is_convex_polygon(square)

        # Closed square
        closed_square = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]]
        assert is_convex_polygon(closed_square)

        # Concave polygon (pac-man or arrowhead)
        concave = [[0.0, 0.0], [2.0, 0.0], [1.0, 1.0], [2.0, 2.0], [0.0, 2.0]]
        assert not is_convex_polygon(concave)

    def test_ccw_and_orientation(self):
        ccw_ring = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.0, 0.0]])
        assert check_counter_clockwise(ccw_ring)
        assert get_polygon_orientation(ccw_ring) == "counter-clockwise"

        cw_ring = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]])
        assert not check_counter_clockwise(cw_ring)
        assert get_polygon_orientation(cw_ring) == "clockwise"

    def test_clean_geometry(self):
        # Self-intersecting bowtie
        wkt = "POLYGON((0 0, 2 2, 2 0, 0 2, 0 0))"
        geom = ogr.CreateGeometryFromWkt(wkt)
        cleaned = clean_geometry(geom)
        assert cleaned is not None
        assert cleaned.IsValid()


class TestRangeAndBoxes:
    def test_longitude_conversions(self):
        assert convert_360_to_180(270.0) == -90.0
        assert convert_360_to_180(0.0) == 0.0
        assert convert_180_to_360(-90.0) == 270.0
        assert convert_180_to_360(0.0) == 0.0

        # Arrays
        arr_360 = np.array([0.0, 90.0, 270.0])
        np.testing.assert_allclose(convert_360_to_180(arr_360), [0.0, 90.0, -90.0])
        arr_180 = np.array([-90.0, 0.0, 90.0])
        np.testing.assert_allclose(convert_180_to_360(arr_180), [270.0, 0.0, 90.0])

    def test_create_box_and_nearest_res(self):
        box, coords = create_box_from_longitude_latitude(-120.0, 45.0, 1.0, 1.0)
        assert len(box) == 4
        assert coords.shape == (4, 2)

        res = find_nearest_resolution(50000.0, 45.0)
        assert res > 0.0


class TestSimplificationAndMesh:
    def test_douglas_peucker(self):
        polyline = [(0.0, 0.0), (0.1, 0.001), (0.2, 0.0), (0.3, 0.0), (0.4, 0.0)]
        simplified = douglas_peucker_geodetic(polyline, tolerance=2000.0)
        assert len(simplified) < len(polyline)

    def test_visvalingam_whyatt(self):
        polyline = [(0.0, 0.0), (0.1, 0.001), (0.2, 0.0), (0.3, 0.0), (0.4, 0.0)]
        simplified = visvalingam_whyatt_geodetic(polyline, tolerance=2000.0)
        assert len(simplified) < len(polyline)

    def test_unique_vertices_and_connectivity(self):
        cell_x = [[0.0, 1.0, 0.5], [1.0, 0.5, 1.5]]
        cell_y = [[0.0, 0.0, 1.0], [0.0, 1.0, 1.0]]
        xv, yv, conn, lookup = extract_unique_vertices_and_connectivity(cell_x, cell_y)
        assert len(xv) == 4
        assert conn.shape == (2, 3)

    def test_output_geometry_type(self):
        assert get_output_geometry_type(ogr.wkbMultiPolygon) == ogr.wkbPolygon
        assert get_output_geometry_type(ogr.wkbLineString25D) == ogr.wkbLineString
