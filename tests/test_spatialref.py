"""Unit tests for pyearth.gis.spatialref package."""

import unittest
import numpy as np
from osgeo import osr

from pyearth.gis.spatialref import (
    get_utm_zone,
    get_utm_epsg_code,
    get_utm_spatial_reference_wkt,
    degree_to_meter,
    meter_to_degree,
    is_wgs84_projection,
    reproject_coordinates,
    reproject_coordinates_batch,
)


class TestUtmUtility(unittest.TestCase):
    def test_get_utm_zone(self):
        # Prime meridian, equator -> Zone 31
        self.assertEqual(get_utm_zone(0.0, 0.0), 31)
        # San Francisco (-122.4, 37.7) -> Zone 10
        self.assertEqual(get_utm_zone(-122.4194, 37.7749), 10)
        # Tokyo (139.7, 35.7) -> Zone 54
        self.assertEqual(get_utm_zone(139.6917, 35.6895), 54)
        # Boundaries: -180 -> Zone 1, 180 -> Zone 60
        self.assertEqual(get_utm_zone(-180.0, 0.0), 1)
        self.assertEqual(get_utm_zone(180.0, 0.0), 60)

    def test_get_utm_zone_validation(self):
        with self.assertRaises(ValueError):
            get_utm_zone(-181.0, 0.0)
        with self.assertRaises(ValueError):
            get_utm_zone(181.0, 0.0)
        with self.assertRaises(ValueError):
            get_utm_zone(0.0, -85.0)
        with self.assertRaises(ValueError):
            get_utm_zone(0.0, 85.0)

    def test_get_utm_epsg_code(self):
        # Northern hemisphere
        self.assertEqual(get_utm_epsg_code(0.0, 51.5), 32631)
        self.assertEqual(get_utm_epsg_code(-122.4, 37.7), 32610)
        # Southern hemisphere
        self.assertEqual(get_utm_epsg_code(-58.38, -34.6), 32721)
        self.assertEqual(get_utm_epsg_code(18.4, -33.9), 32734)

    def test_get_utm_spatial_reference_wkt(self):
        wkt = get_utm_spatial_reference_wkt(-122.4194, 37.7749)
        self.assertIsInstance(wkt, str)
        self.assertTrue("UTM" in wkt or "32610" in wkt)


class TestDegreeMeterConversion(unittest.TestCase):
    def test_degree_to_meter_scalar(self):
        # At equator, 1 degree is roughly 111.3 km
        m_equator = degree_to_meter(1.0, 0.0)
        self.assertAlmostEqual(m_equator, 111319.49, delta=1.0)

        # At 60 degrees, 1 degree longitude is half the equatorial value
        m_60 = degree_to_meter(1.0, 60.0)
        self.assertAlmostEqual(m_60, m_equator / 2.0, delta=1.0)

    def test_meter_to_degree_scalar(self):
        m_equator = 111319.49079327358
        deg = meter_to_degree(m_equator, 0.0)
        self.assertAlmostEqual(deg, 1.0, places=5)

    def test_conversion_roundtrip(self):
        original_deg = 0.25
        lat = 35.5
        m = degree_to_meter(original_deg, lat)
        recovered_deg = meter_to_degree(m, lat)
        self.assertAlmostEqual(original_deg, recovered_deg, places=7)

    def test_vectorized_conversion(self):
        degs = np.array([0.5, 1.0, 2.0])
        lats = np.array([0.0, 30.0, 60.0])
        meters = degree_to_meter(degs, lats)
        self.assertIsInstance(meters, np.ndarray)
        self.assertEqual(len(meters), 3)

        recovered = meter_to_degree(meters, lats)
        np.testing.assert_allclose(recovered, degs, rtol=1e-6)

    def test_polar_exception(self):
        with self.assertRaises(ValueError):
            meter_to_degree(100.0, 90.0)
        with self.assertRaises(ValueError):
            meter_to_degree(100.0, -90.0)

    def test_legacy_kwargs(self):
        m = degree_to_meter(dResolution_degree_in=1.0, dLatitude_mean_in=0.0)
        self.assertAlmostEqual(m, 111319.49, delta=1.0)
        d = meter_to_degree(dResolution_meter_in=m, dLatitude_mean_in=0.0)
        self.assertAlmostEqual(d, 1.0, places=5)


class TestIsWgs84Projection(unittest.TestCase):
    def test_wgs84_variants(self):
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(4326)
        wkt = sr.ExportToWkt()

        self.assertTrue(is_wgs84_projection("EPSG:4326"))
        self.assertTrue(is_wgs84_projection("4326"))
        self.assertTrue(is_wgs84_projection(wkt))
        self.assertTrue(is_wgs84_projection(sr))

    def test_non_wgs84(self):
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(3857)
        self.assertFalse(is_wgs84_projection(sr.ExportToWkt()))
        self.assertFalse(is_wgs84_projection(sr))
        self.assertFalse(is_wgs84_projection("EPSG:3857"))

    def test_invalid_inputs(self):
        self.assertFalse(is_wgs84_projection(None))
        self.assertFalse(is_wgs84_projection(""))
        self.assertFalse(is_wgs84_projection("not_a_valid_projection_definition"))


class TestReprojectCoordinates(unittest.TestCase):
    def test_reproject_single(self):
        # Reproject (0.0, 0.0) in WGS84 to UTM Zone 31N
        x_utm, y_utm = reproject_coordinates(0.0, 0.0, "EPSG:4326", "EPSG:32631")
        self.assertAlmostEqual(x_utm, 166021.44, delta=1.0)
        self.assertAlmostEqual(y_utm, 0.0, delta=1.0)

        # Reproject back to WGS84 (using default target)
        x_wgs, y_wgs = reproject_coordinates(x_utm, y_utm, "EPSG:32631")
        self.assertAlmostEqual(x_wgs, 0.0, delta=1e-5)
        self.assertAlmostEqual(y_wgs, 0.0, delta=1e-5)

    def test_reproject_batch(self):
        xs = [0.0, 1.0, 2.0]
        ys = [0.0, 1.0, 2.0]
        xs_utm, ys_utm = reproject_coordinates_batch(xs, ys, "EPSG:4326", "EPSG:32631")
        self.assertEqual(len(xs_utm), 3)
        self.assertEqual(len(ys_utm), 3)

        # Batch back to WGS84
        xs_wgs, ys_wgs = reproject_coordinates_batch(xs_utm, ys_utm, "EPSG:32631")
        for orig_x, rec_x in zip(xs, xs_wgs):
            self.assertAlmostEqual(orig_x, rec_x, delta=1e-5)

    def test_reproject_validation(self):
        with self.assertRaises(ValueError):
            reproject_coordinates(None, 0.0, "EPSG:4326")
        with self.assertRaises(ValueError):
            reproject_coordinates_batch([], [], "EPSG:4326")
        with self.assertRaises(ValueError):
            reproject_coordinates_batch([1.0], [1.0, 2.0], "EPSG:4326")


if __name__ == "__main__":
    unittest.main()
