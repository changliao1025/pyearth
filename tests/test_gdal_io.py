"""Unit tests for GDAL raster and vector I/O utilities."""

import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest
from osgeo import gdal, ogr, osr

from pyearth.gis.gdal.read import (
    gdal_get_raster_extent,
    gdal_get_raster_spatial_reference,
    gdal_get_raster_spatial_reference_wkt,
    gdal_get_vector_boundary,
    gdal_get_vector_extent,
    gdal_get_vector_spatial_reference,
    gdal_read_ascii_file,
    gdal_read_envi_file,
    gdal_read_envi_file_multiple_band,
    gdal_read_geotiff_file,
    gdal_read_geotiff_file_multiple_band,
    gdal_read_vector_file,
)
from pyearth.gis.gdal.write import (
    export_point_as_polygon_file,
    export_point_to_vector_file,
    gdal_write_envi_file,
    gdal_write_envi_file_multiple_band,
    gdal_write_geotiff_file,
    gdal_write_geotiff_file_multiple_band,
    gdal_write_wkt_to_vector_file,
)


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test artifacts and clean up after."""
    d = tempfile.mkdtemp()
    yield Path(d)
    shutil.rmtree(d, ignore_errors=True)


class DummyPoint:
    """Mock point object with coordinates."""

    def __init__(self, lon, lat, x=None, y=None):
        self.dLongitude_degree = lon
        self.dLatitude_degree = lat
        self.dx = x if x is not None else lon
        self.dy = y if y is not None else lat


class TestGDALRasterIO:
    """Tests for raster reading and writing."""

    def test_geotiff_single_band_roundtrip(self, temp_dir):
        out_tif = temp_dir / "test_single.tif"
        data = np.arange(100, dtype=np.float32).reshape((10, 10))

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        gdal_write_geotiff_file(
            filename=out_tif,
            data=data,
            pixel_width=0.1,
            pixel_height=-0.1,
            origin_x=-120.0,
            origin_y=45.0,
            missing_value=-9999.0,
            projection=srs,
        )

        assert out_tif.exists()

        result = gdal_read_geotiff_file(out_tif)
        read_data = result["data"]
        geotransform = result["geotransform"]
        missing_val = result["missing_value"]

        np.testing.assert_allclose(read_data, data)
        assert missing_val == -9999.0
        assert pytest.approx(geotransform[0]) == -120.0
        assert pytest.approx(geotransform[1]) == 0.1
        assert pytest.approx(geotransform[3]) == 45.0
        assert pytest.approx(geotransform[5]) == -0.1

        # Test raster spatial reference
        srs_read = gdal_get_raster_spatial_reference(out_tif)
        assert srs_read is not None
        wkt_read = gdal_get_raster_spatial_reference_wkt(out_tif)
        assert len(wkt_read) > 0

        # Test raster extent
        extent = gdal_get_raster_extent(out_tif)
        assert pytest.approx(extent[0]) == -120.0
        assert pytest.approx(extent[1]) == -119.0
        assert pytest.approx(extent[2]) == 44.0
        assert pytest.approx(extent[3]) == 45.0

    def test_geotiff_multiband_roundtrip(self, temp_dir):
        out_tif = temp_dir / "test_multi.tif"
        data = np.stack(
            [
                np.ones((10, 10), dtype=np.int32) * 1,
                np.ones((10, 10), dtype=np.int32) * 2,
                np.ones((10, 10), dtype=np.int32) * 3,
            ],
            axis=0,
        )

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(3857)

        gdal_write_geotiff_file_multiple_band(
            filename=out_tif,
            data=data,
            pixel_width=10.0,
            pixel_height=10.0,
            origin_x=0.0,
            origin_y=0.0,
            projection=srs,
        )

        assert out_tif.exists()

        result = gdal_read_geotiff_file_multiple_band(out_tif)
        read_data = result["data"]
        assert read_data.shape == (3, 10, 10)
        np.testing.assert_array_equal(read_data, data)

    def test_envi_single_and_multiband_roundtrip(self, temp_dir):
        out_bin = temp_dir / "test_envi.dat"
        data = np.arange(25, dtype=np.float32).reshape((5, 5))

        gdal_write_envi_file(
            filename=out_bin,
            data=data,
            pixel_width=1.0,
            pixel_height=1.0,
            origin_x=10.0,
            origin_y=20.0,
            missing_value=-999.0,
        )

        result = gdal_read_envi_file(out_bin)
        np.testing.assert_allclose(result["data"], data)

        # Multi-band ENVI
        out_multi_bin = temp_dir / "test_envi_multi.dat"
        data_multi = np.ones((2, 5, 5), dtype=np.float32) * 5.0
        gdal_write_envi_file_multiple_band(
            filename=out_multi_bin,
            data=data_multi,
            pixel_width=1.0,
            pixel_height=1.0,
            origin_x=10.0,
            origin_y=20.0,
        )

        result_multi = gdal_read_envi_file_multiple_band(out_multi_bin)
        assert result_multi["data"].shape == (2, 5, 5)
        np.testing.assert_allclose(result_multi["data"], data_multi)

    def test_ascii_read(self, temp_dir):
        asc_path = temp_dir / "test.asc"
        content = (
            "ncols         4\n"
            "nrows         4\n"
            "xllcorner     0.0\n"
            "yllcorner     0.0\n"
            "cellsize      1.0\n"
            "NODATA_value  -9999\n"
            "1 2 3 4\n"
            "5 6 7 8\n"
            "9 10 11 12\n"
            "13 14 15 16\n"
        )
        asc_path.write_text(content)

        result = gdal_read_ascii_file(asc_path)
        assert result["data"].shape == (4, 4)
        assert result["data"][0, 0] == 1.0
        assert result["missing_value"] == -9999.0


class TestGDALVectorIO:
    """Tests for vector reading and writing."""

    def test_wkt_write_and_read(self, temp_dir):
        geojson_path = temp_dir / "polygon.geojson"
        wkt = "POLYGON ((-120 40, -119 40, -119 41, -120 41, -120 40))"

        gdal_write_wkt_to_vector_file(
            wkt=wkt,
            sFilename_out=geojson_path,
            spatial_reference=4326,
            overwrite=True,
        )

        assert geojson_path.exists()

        result = gdal_read_vector_file(geojson_path)
        assert result["feature_count"] == 1
        assert len(result["geometries"]) == 1

        # Check vector extent and boundary
        extent = gdal_get_vector_extent(geojson_path)
        assert pytest.approx(extent[0]) == -120.0
        assert pytest.approx(extent[1]) == -119.0
        assert pytest.approx(extent[2]) == 40.0
        assert pytest.approx(extent[3]) == 41.0

        boundary_wkt, (min_x, max_x, min_y, max_y) = gdal_get_vector_boundary(
            geojson_path
        )
        assert pytest.approx(min_x) == -120.0
        assert pytest.approx(max_x) == -119.0

        srs = gdal_get_vector_spatial_reference(geojson_path)
        assert srs is not None

    def test_export_point_objects_and_tuples(self, temp_dir):
        geojson_path = temp_dir / "points.geojson"
        points = [
            DummyPoint(-122.0, 37.0),
            DummyPoint(-121.0, 38.0),
            {"dLongitude_degree": -120.0, "dLatitude_degree": 39.0},
            (-119.0, 40.0),
        ]
        attrs = [10, 20, 30, 40]

        export_point_to_vector_file(
            aPoint_in=points,
            sFilename_vector_in=geojson_path,
            iFlag_projected_in=0,
            pSpatial_reference_in=4326,
            aAttribute_data=attrs,
            overwrite=True,
        )

        assert geojson_path.exists()

        result = gdal_read_vector_file(geojson_path)
        assert result["feature_count"] == 4
        assert len(result["attributes"]) == 4
        assert result["attributes"][0]["pointid"] == 1
        assert result["attributes"][0]["connectivity"] == 10
        assert result["attributes"][3]["pointid"] == 4
        assert result["attributes"][3]["connectivity"] == 40

    def test_export_point_as_polygon(self, temp_dir):
        shp_path = temp_dir / "test_poly.shp"
        points = [
            DummyPoint(0.0, 0.0),
            DummyPoint(1.0, 0.0),
            DummyPoint(1.0, 1.0),
            DummyPoint(0.0, 1.0),
        ]

        export_point_as_polygon_file(
            aPoint_in=points,
            sFilename_out=shp_path,
            pSpatial_reference_in=4326,
            overwrite=True,
        )

        assert shp_path.exists()

        result = gdal_read_vector_file(shp_path)
        assert result["feature_count"] == 1
        assert "area" in result["attributes"][0]
        assert result["attributes"][0]["area"] > 0
