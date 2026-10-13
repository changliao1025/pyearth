import numpy as np
import pytest
from osgeo import ogr

from pyearth.gis.idl_handler import (
    IdlHandler,
    IdlStrategy,
    IdlResult,
    check_cross_international_date_line_geometry,
    check_cross_international_date_line_polygon,
    split_international_date_line_polygon_coordinates,
    unwrap_longitudes,
)


def test_idl_handler_detect_crossing():
    # Crossing polygon
    coords = np.array([
        [170.0, -10.0],
        [170.0,  10.0],
        [-170.0, 10.0],
        [-170.0,-10.0],
        [170.0, -10.0],
    ])
    handler = IdlHandler(strategy=IdlStrategy.CHECK_ONLY)
    result = handler.detect(coords)
    assert result.crosses_idl is True
    assert result.sub_polygons is None


def test_idl_handler_detect_non_crossing():
    coords = np.array([
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
        [0.0, 0.0],
    ])
    handler = IdlHandler(strategy=IdlStrategy.CHECK_ONLY)
    result = handler.detect(coords)
    assert result.crosses_idl is False


def test_idl_handler_split_strategy():
    coords = np.array([
        [170.0, -10.0],
        [170.0,  10.0],
        [-170.0, 10.0],
        [-170.0,-10.0],
        [170.0, -10.0],
    ])
    handler = IdlHandler(strategy=IdlStrategy.SPLIT)
    result = handler.detect(coords)
    assert result.crosses_idl is True
    assert result.sub_polygons is not None
    assert len(result.sub_polygons) == 2

    eastern, western = result.sub_polygons
    assert np.all(eastern[:, 0] > 0)
    assert np.all(western[:, 0] < 0)
    assert np.allclose(eastern[0], eastern[-1])
    assert np.allclose(western[0], western[-1])


def test_idl_handler_split_starting_western():
    # Polygon starting with negative longitude
    coords = np.array([
        [-170.0, -10.0],
        [170.0, -10.0],
        [170.0,  10.0],
        [-170.0, 10.0],
        [-170.0, -10.0],
    ])
    sub_polys = split_international_date_line_polygon_coordinates(coords)
    assert len(sub_polys) == 2
    eastern, western = sub_polys
    assert np.all(eastern[:, 0] > 0)
    assert np.all(western[:, 0] < 0)
    assert np.allclose(eastern[0], eastern[-1])
    assert np.allclose(western[0], western[-1])


def test_idl_handler_unwrap():
    coords = np.array([
        [170.0, 10.0],
        [-170.0, 10.0],
        [-170.0, -10.0],
        [170.0, -10.0],
        [170.0, 10.0],
    ])
    handler = IdlHandler()
    unwrapped = handler.unwrap(coords)
    diffs = np.abs(np.diff(unwrapped[:, 0]))
    assert np.all(diffs < 180.0)


def test_idl_handler_check_ccw():
    # Counter-clockwise polygon across IDL
    coords_ccw = np.array([
        [170.0, -10.0],
        [-170.0, -10.0],
        [-170.0, 10.0],
        [170.0, 10.0],
        [170.0, -10.0],
    ])
    # Shoelace on unwrapped should determine orientation
    assert IdlHandler.check_counter_clockwise(coords_ccw) is True
    assert IdlHandler.check_counter_clockwise(coords_ccw[::-1]) is False


def test_idl_handler_fix_ogr_geometry():
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(170.0, -10.0)
    ring.AddPoint(170.0, 10.0)
    ring.AddPoint(-170.0, 10.0)
    ring.AddPoint(-170.0, -10.0)
    ring.AddPoint(170.0, -10.0)
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)

    assert check_cross_international_date_line_geometry(poly) is True

    # Test SPLIT
    split_handler = IdlHandler(strategy=IdlStrategy.SPLIT)
    fixed_geoms = split_handler.fix(poly)
    assert len(fixed_geoms) == 2
    for g in fixed_geoms:
        assert g.IsValid()

    # Test DROP
    drop_handler = IdlHandler(strategy=IdlStrategy.DROP)
    assert drop_handler.fix(poly) == []

    # Test WRAP_TO_360
    wrap_handler = IdlHandler(strategy=IdlStrategy.WRAP_TO_360)
    wrapped_geoms = wrap_handler.fix(poly)
    assert len(wrapped_geoms) == 1
    assert wrapped_geoms[0].IsValid()
