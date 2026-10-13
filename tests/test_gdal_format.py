"""Unit tests for pyearth.gis.envi and pyearth.gis.gdal format utilities."""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from osgeo import gdal, gdalconst, ogr

from pyearth.gis.envi import envi_write_header
from pyearth.gis.gdal import (
    gdal_check_file_type,
    gdal_check_raster_valid,
    gdal_raster_format_support,
    gdal_to_numpy_datatype,
    gdal_validate_polygon_file,
    gdal_vector_format_support,
    get_available_raster_formats,
    get_available_vector_formats,
    get_extension_from_raster_format,
    get_extension_from_vector_format,
    get_raster_driver_from_filename,
    get_raster_driver_from_format,
    get_raster_format_from_extension,
    get_raster_format_from_filename,
    get_vector_driver_from_filename,
    get_vector_driver_from_format,
    get_vector_format_from_extension,
    get_vector_format_from_filename,
    numpy_dtype_to_gdal,
    numpy_to_gdal_type,
)


class TestEnviWriteHeader(unittest.TestCase):
    def test_write_valid_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            hdr_path = Path(tmpdir) / "test.hdr"
            meta = {
                "sFilename": "test.dat",
                "ncolumn": 50,
                "nrow": 60,
                "nband": 2,
                "offset": 0,
                "data_type": 4,
                "byte_order": 0,
                "missing_value": -9999,
                "ULlon": -120.0,
                "ULlat": 38.0,
                "pixelSize": 0.05,
                "interleave": "bsq",
            }
            envi_write_header(hdr_path, meta)
            self.assertTrue(hdr_path.exists())

            content = hdr_path.read_text(encoding="utf-8")
            self.assertTrue(content.startswith("ENVI\n"))
            self.assertIn("samples = 50", content)
            self.assertIn("lines = 60", content)
            self.assertIn("bands = 2", content)
            self.assertIn("interleave = bsq", content)
            self.assertIn("WGS-84", content)

    def test_missing_keys(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            hdr_path = Path(tmpdir) / "test.hdr"
            with self.assertRaises(ValueError):
                envi_write_header(hdr_path, {"sFilename": "incomplete.dat"})


class TestGdalRasterFormatSupport(unittest.TestCase):
    def test_format_mappings(self):
        self.assertEqual(get_raster_format_from_extension(".tif"), "GTiff")
        self.assertEqual(get_raster_format_from_extension("tiff"), "GTiff")
        self.assertEqual(get_raster_format_from_extension(".png"), "PNG")
        self.assertEqual(get_raster_format_from_extension(".nc"), "netCDF")

        with self.assertRaises(ValueError):
            get_raster_format_from_extension(".unknown_raster_ext")

    def test_format_from_filename(self):
        self.assertEqual(get_raster_format_from_filename("path/to/image.tif"), "GTiff")
        self.assertEqual(get_raster_format_from_filename(Path("image.png")), "PNG")

        with self.assertRaises(ValueError):
            get_raster_format_from_filename("no_extension")

    def test_reverse_lookup(self):
        self.assertEqual(get_extension_from_raster_format("GTiff"), ".tif")
        self.assertEqual(get_extension_from_raster_format("PNG"), ".png")

    def test_driver_lookup(self):
        drv = get_raster_driver_from_format("GTiff")
        self.assertIsNotNone(drv)
        self.assertEqual(drv.GetName(), "GTiff")

        drv_file = get_raster_driver_from_filename("test.png")
        self.assertEqual(drv_file.GetName(), "PNG")

    def test_available_formats(self):
        avail = get_available_raster_formats()
        self.assertIsInstance(avail, dict)
        self.assertIn(".tif", avail)
        self.assertEqual(gdal_raster_format_support(), avail)


class TestGdalVectorFormatSupport(unittest.TestCase):
    def test_format_mappings(self):
        self.assertEqual(get_vector_format_from_extension(".shp"), "ESRI Shapefile")
        self.assertEqual(get_vector_format_from_extension(".geojson"), "GeoJSON")
        self.assertEqual(get_vector_format_from_extension(".gpkg"), "GPKG")

        with self.assertRaises(ValueError):
            get_vector_format_from_extension(".unknown_vector_ext")

    def test_format_from_filename(self):
        self.assertEqual(get_vector_format_from_filename("vector.geojson"), "GeoJSON")
        self.assertEqual(get_vector_format_from_filename(Path("data.shp")), "ESRI Shapefile")

        with self.assertRaises(ValueError):
            get_vector_format_from_filename("no_extension")

    def test_reverse_lookup(self):
        self.assertEqual(get_extension_from_vector_format("GeoJSON"), ".geojson")
        self.assertEqual(get_extension_from_vector_format("ESRI Shapefile"), ".shp")

    def test_driver_lookup(self):
        drv = get_vector_driver_from_format("GeoJSON")
        self.assertIsNotNone(drv)
        self.assertEqual(drv.GetName(), "GeoJSON")

        drv_file = get_vector_driver_from_filename("roads.shp")
        self.assertEqual(drv_file.GetName(), "ESRI Shapefile")

    def test_available_formats(self):
        avail = get_available_vector_formats()
        self.assertIsInstance(avail, dict)
        self.assertIn(".geojson", avail)
        self.assertEqual(gdal_vector_format_support(), avail)


class TestGdalTypeConversions(unittest.TestCase):
    def test_gdal_to_numpy(self):
        self.assertEqual(gdal_to_numpy_datatype(gdalconst.GDT_Byte), np.uint8)
        self.assertEqual(gdal_to_numpy_datatype(gdalconst.GDT_Int16), np.int16)
        self.assertEqual(gdal_to_numpy_datatype(gdalconst.GDT_Int32), np.int32)
        self.assertEqual(gdal_to_numpy_datatype(gdalconst.GDT_Float32), np.float32)
        self.assertEqual(gdal_to_numpy_datatype(gdalconst.GDT_Float64), np.float64)

        with self.assertRaises(ValueError):
            gdal_to_numpy_datatype(999999)

    def test_numpy_to_gdal(self):
        self.assertEqual(numpy_dtype_to_gdal(np.uint8), gdalconst.GDT_Byte)
        self.assertEqual(numpy_dtype_to_gdal(np.int16), gdalconst.GDT_Int16)
        self.assertEqual(numpy_dtype_to_gdal(np.float32), gdalconst.GDT_Float32)
        self.assertEqual(numpy_dtype_to_gdal(np.float64), gdalconst.GDT_Float64)
        self.assertEqual(numpy_dtype_to_gdal("int32"), gdalconst.GDT_Int32)

    def test_numpy_to_gdal_value(self):
        self.assertEqual(numpy_to_gdal_type(np.int32(42)), 42)
        self.assertEqual(numpy_to_gdal_type(np.float64(3.14)), 3.14)
        self.assertEqual(numpy_to_gdal_type(np.nan, target_type=ogr.OFTReal), 0.0)
        self.assertEqual(numpy_to_gdal_type(None, target_type=ogr.OFTString), "")


class TestDatasetInspectionAndValidation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.raster_path = Path(self.tmpdir.name) / "test.tif"
        self.vector_path = Path(self.tmpdir.name) / "test.geojson"

        # Create dummy GeoTIFF
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(str(self.raster_path), 10, 10, 1, gdalconst.GDT_Byte)
        band = ds.GetRasterBand(1)
        band.WriteArray(np.zeros((10, 10), dtype=np.uint8))
        ds.FlushCache()
        ds = None

        # Create dummy GeoJSON with a valid polygon
        ogr_driver = ogr.GetDriverByName("GeoJSON")
        vds = ogr_driver.CreateDataSource(str(self.vector_path))
        layer = vds.CreateLayer("polys", geom_type=ogr.wkbPolygon)
        feat_defn = layer.GetLayerDefn()
        feat = ogr.Feature(feat_defn)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(0.0, 0.0)
        ring.AddPoint(1.0, 0.0)
        ring.AddPoint(1.0, 1.0)
        ring.AddPoint(0.0, 1.0)
        ring.AddPoint(0.0, 0.0)
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feat.SetGeometry(poly)
        layer.CreateFeature(feat)
        vds = None

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_gdal_check_file_type(self):
        self.assertEqual(gdal_check_file_type(self.raster_path), "raster")
        self.assertEqual(gdal_check_file_type(self.vector_path), "vector")
        self.assertEqual(gdal_check_file_type("/non/existent.dat"), "unknown")

    def test_gdal_check_raster_valid(self):
        self.assertTrue(gdal_check_raster_valid(self.raster_path))
        self.assertFalse(gdal_check_raster_valid(self.vector_path))
        self.assertFalse(gdal_check_raster_valid("/non/existent.tif"))

    def test_gdal_validate_polygon_file(self):
        self.assertTrue(gdal_validate_polygon_file(self.vector_path))
        self.assertFalse(gdal_validate_polygon_file(self.raster_path))
        self.assertFalse(gdal_validate_polygon_file("/non/existent.shp"))


if __name__ == "__main__":
    unittest.main()
