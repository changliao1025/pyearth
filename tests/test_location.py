"""Unit tests for GIS location, coordinate conversion, and grid indexing utilities."""

import numpy as np
import pytest
from osgeo import ogr

from pyearth.gis.location import (
    convert_longitude_latitude_to_sphere_3d,
    convert_sphere_3d_to_longitude_latitude,
    find_index_by_longitude_latitude,
    get_geometry_coordinates,
    get_hydrosheds_continent_from_extent,
    get_linearring_coords,
    get_linestring_coords,
    get_multipolygon_exterior_coords,
    get_point_coords,
    get_polygon_exterior_coords,
    google_meters_per_pixel,
    pixel_to_world,
    world_to_pixel,
)


class TestGoogleMetersPerPixel:
    def test_equatorial_resolutions(self):
        res0 = google_meters_per_pixel(0)
        assert pytest.approx(res0, rel=1e-3) == 156543.03

        res10 = google_meters_per_pixel(10)
        assert pytest.approx(res10, rel=1e-3) == 152.87

    def test_latitude_scaling(self):
        res_eq = google_meters_per_pixel(10)
        res_60 = google_meters_per_pixel(10, latitude=60.0)
        # cos(60 deg) = 0.5
        assert pytest.approx(res_60, rel=1e-3) == res_eq * 0.5

    def test_invalid_zoom_and_latitude(self):
        with pytest.raises(ValueError):
            google_meters_per_pixel(-1)
        with pytest.raises(ValueError):
            google_meters_per_pixel(31)
        with pytest.raises(ValueError):
            google_meters_per_pixel(10, latitude=95.0)


class TestWorldAndPixel:
    def test_north_up_roundtrip(self):
        # origin (-120.0, 45.0), pixel width 0.1, pixel height -0.1
        gt = [-120.0, 0.1, 0.0, 45.0, 0.0, -0.1]

        col, row = world_to_pixel(gt, -119.5, 44.5)
        assert col == 5
        assert row == 5

        # Pixel to world (top-left)
        x, y = pixel_to_world(gt, col, row, center=False)
        assert pytest.approx(x) == -119.5
        assert pytest.approx(y) == 44.5

        # Pixel to world (center)
        xc, yc = pixel_to_world(gt, col, row, center=True)
        assert pytest.approx(xc) == -119.45
        assert pytest.approx(yc) == 44.45

    def test_rotated_affine_roundtrip(self):
        # 45-degree rotation
        c = np.cos(np.radians(45))
        s = np.sin(np.radians(45))
        gt = [100.0, c, -s, 200.0, s, c]

        col, row = world_to_pixel(gt, 105.0, 205.0)
        x_rec, y_rec = pixel_to_world(gt, col, row)
        # Should be within 1 pixel
        assert abs(x_rec - 105.0) < 2.0
        assert abs(y_rec - 205.0) < 2.0


class TestSphere3DConversion:
    def test_cardinal_points(self):
        # North pole
        np_xyz = convert_longitude_latitude_to_sphere_3d(0.0, 90.0)
        np.testing.assert_allclose(np_xyz, [0.0, 0.0, 1.0], atol=1e-10)

        # Prime meridian equator
        pm_xyz = convert_longitude_latitude_to_sphere_3d(0.0, 0.0)
        np.testing.assert_allclose(pm_xyz, [1.0, 0.0, 0.0], atol=1e-10)

        # 90E equator
        e90_xyz = convert_longitude_latitude_to_sphere_3d(90.0, 0.0)
        np.testing.assert_allclose(e90_xyz, [0.0, 1.0, 0.0], atol=1e-10)

    def test_roundtrip(self):
        test_coords = [(45.0, 30.0), (-120.0, -45.0), (179.0, 89.0)]
        for lon, lat in test_coords:
            xyz = convert_longitude_latitude_to_sphere_3d(lon, lat)
            lon_rec, lat_rec = convert_sphere_3d_to_longitude_latitude(*xyz)
            assert pytest.approx(lon_rec, abs=1e-7) == lon
            assert pytest.approx(lat_rec, abs=1e-7) == lat

    def test_custom_radius_and_vectorization(self):
        lons = np.array([0.0, 90.0])
        lats = np.array([0.0, 0.0])
        xyz_arr = convert_longitude_latitude_to_sphere_3d(lons, lats, radius=10.0)
        assert xyz_arr.shape == (2, 3)
        np.testing.assert_allclose(xyz_arr[0], [10.0, 0.0, 0.0], atol=1e-10)
        np.testing.assert_allclose(xyz_arr[1], [0.0, 10.0, 0.0], atol=1e-10)


class TestFindIndexByLonLat:
    def test_center_and_corner_modes(self):
        # Cell-center registration
        row, col = find_index_by_longitude_latitude(
            dLon_in=3.5,
            dLat_in=47.5,
            dLon_min_in=0.0,
            dLat_max_in=50.0,
            dResolution_x_in=1.0,
            dResolution_y_in=1.0,
            iFlag_center_in=True,
        )
        assert row == 2
        assert col == 4

        # Corner registration
        row_c, col_c = find_index_by_longitude_latitude(
            dLon_in=3.5,
            dLat_in=47.5,
            dLon_min_in=-0.5,
            dLat_max_in=50.5,
            dResolution_x_in=1.0,
            dResolution_y_in=1.0,
            iFlag_center_in=False,
        )
        assert row_c == 2
        assert col_c == 4

    def test_out_of_bounds_errors(self):
        with pytest.raises(ValueError):
            find_index_by_longitude_latitude(
                dLon_in=10.0,
                dLat_in=60.0,  # North of 50
                dLon_min_in=0.0,
                dLat_max_in=50.0,
                dResolution_x_in=1.0,
                dResolution_y_in=1.0,
            )


class TestHydroSHEDSContinent:
    def test_north_america(self):
        # Texas
        extent = (-106.0, -93.5, 25.8, 36.5)
        cont = get_hydrosheds_continent_from_extent(extent)
        assert cont == "na"

    def test_europe(self):
        extent = (-10.0, 30.0, 35.0, 60.0)
        cont = get_hydrosheds_continent_from_extent(extent)
        assert cont == "eu"

    def test_return_all(self):
        extent = (30.0, 80.0, 40.0, 60.0)
        conts = get_hydrosheds_continent_from_extent(extent, return_all=True)
        assert isinstance(conts, list)
        assert "eu" in conts or "as" in conts


class TestGeometryCoordinates:
    def test_point_and_linestring(self):
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(10.0, 20.0)
        coords = get_geometry_coordinates(pt)
        assert coords.shape == (1, 2)
        assert coords[0, 0] == 10.0
        assert coords[0, 1] == 20.0

        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(0.0, 0.0)
        line.AddPoint(1.0, 1.0)
        line_coords = get_geometry_coordinates(line)
        assert line_coords.shape == (2, 2)

    def test_polygon_ccw_enforcement(self):
        # Clockwise square
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0.0, 0.0)
        ring.AddPoint(0.0, 1.0)
        ring.AddPoint(1.0, 1.0)
        ring.AddPoint(1.0, 0.0)
        ring.AddPoint(0.0, 0.0)

        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)

        coords_original = get_geometry_coordinates(poly, enforce_ccw=False)
        assert coords_original.shape == (5, 2)

        coords_ccw = get_geometry_coordinates(poly, enforce_ccw=True)
        assert coords_ccw.shape == (5, 2)
        # When reversed, second point should be (1.0, 0.0)
        assert coords_ccw[1, 0] == 1.0
        assert coords_ccw[1, 1] == 0.0
