"""Unit tests for remove_internal_polygon."""

import json
import os
import tempfile

import pytest
from osgeo import ogr

from pyearth.toolbox.management.vector import remove_internal_polygon


def test_remove_internal_polygon():
    geojson_data = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {"id": 1, "name": "poly_with_hole"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [0.0, 0.0],
                            [10.0, 0.0],
                            [10.0, 10.0],
                            [0.0, 10.0],
                            [0.0, 0.0],
                        ],
                        [[2.0, 2.0], [8.0, 2.0], [8.0, 8.0], [2.0, 8.0], [2.0, 2.0]],
                    ],
                },
            },
            {
                "type": "Feature",
                "properties": {"id": 2, "name": "poly_without_hole"},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [20.0, 20.0],
                            [30.0, 20.0],
                            [30.0, 30.0],
                            [20.0, 30.0],
                            [20.0, 20.0],
                        ]
                    ],
                },
            },
        ],
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = os.path.join(tmpdir, "in.geojson")
        out_file = os.path.join(tmpdir, "out.geojson")
        with open(in_file, "w") as f:
            json.dump(geojson_data, f)

        remove_internal_polygon(in_file, out_file, verbose=False)

        ds = ogr.Open(out_file)
        assert ds is not None
        layer = ds.GetLayer()
        assert layer.GetFeatureCount() == 2

        f1 = layer.GetNextFeature()
        g1 = f1.GetGeometryRef()
        assert g1.GetGeometryCount() == 1  # Hole removed!
        assert f1.GetField("name") == "poly_with_hole"

        f2 = layer.GetNextFeature()
        g2 = f2.GetGeometryRef()
        assert g2.GetGeometryCount() == 1
        assert f2.GetField("name") == "poly_without_hole"
