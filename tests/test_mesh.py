import numpy as np
import pytest

from pyearth.toolbox.mesh import (
    pycircle,
    pyline,
    pynvector,
    pypoint,
    pypolygon,
    pypolyline,
)
from pyearth.toolbox.mesh.algorithm import (
    nvector_to_point,
    point_to_nvector,
    slerp,
    split_line_by_length,
)


class TestPyPoint:
    def test_init_and_properties(self):
        pt = pypoint({"dLongitude_degree": 120.0, "dLatitude_degree": 30.0})
        assert pt.dLongitude_degree == 120.0
        assert pt.dLatitude_degree == 30.0
        assert np.isclose(pt.dLongitude_radian, np.radians(120.0))
        assert np.isclose(pt.dLatitude_radian, np.radians(30.0))

    def test_wkt_and_json(self):
        pt = pypoint({"dLongitude_degree": 0.0, "dLatitude_degree": 0.0})
        wkt = pt.towkt()
        assert "POINT (0.0 0.0)" in wkt
        js = pt.tojson()
        assert "dLongitude_degree" in js

    def test_calculate_xyz(self):
        pt = pypoint({"dLongitude_degree": 0.0, "dLatitude_degree": 0.0})
        x, y, z = pt.calculate_xyz()
        assert np.isclose(x, 6371000.0)
        assert np.isclose(y, 0.0)
        assert np.isclose(z, 0.0)


class TestPyNvector:
    def test_init_and_length(self):
        nv = pynvector({"x": 1.0, "y": 0.0, "z": 0.0})
        assert np.isclose(nv.length(), 1.0)

    def test_operations(self):
        nv1 = pynvector({"x": 1.0, "y": 0.0, "z": 0.0})
        nv2 = pynvector({"x": 0.0, "y": 1.0, "z": 0.0})
        assert np.isclose(nv1.dot(nv2), 0.0)

        nv3 = (nv1 + nv2).normalize()
        assert np.isclose(nv3.length(), 1.0)
        assert np.isclose(nv3.dX, nv3.dY)

    def test_high_precision(self):
        nv = pynvector({"x": 0.5, "y": 0.5, "z": 0.7071}, use_high_precision=True)
        assert nv.use_high_precision is True
        assert np.isclose(nv.length(), 1.0)


class TestNvectorPointAdapter:
    def test_round_trip(self):
        pt = pypoint({"dLongitude_degree": 45.0, "dLatitude_degree": -30.0})
        nv = point_to_nvector(pt)
        pt_back = nvector_to_point(nv)
        assert np.isclose(pt_back.dLongitude_degree, 45.0)
        assert np.isclose(pt_back.dLatitude_degree, -30.0)


class TestPyLine:
    def test_line_init(self):
        p1 = pypoint({"dLongitude_degree": 0.0, "dLatitude_degree": 0.0})
        p2 = pypoint({"dLongitude_degree": 1.0, "dLatitude_degree": 0.0})
        line = pyline(p1, p2)
        assert line.pPoint_start == p1
        assert line.pPoint_end == p2
        assert line.calculate_length() > 0.0
        js = line.tojson()
        assert "pPoint_start" in js

    def test_slerp_and_split(self):
        p1 = pypoint({"dLongitude_degree": 0.0, "dLatitude_degree": 0.0})
        p2 = pypoint({"dLongitude_degree": 0.0, "dLatitude_degree": 10.0})
        line = pyline(p1, p2)
        # Split every 200 km (200000 m)
        segments = split_line_by_length(line, 200000.0)
        assert len(segments) >= 3
        assert np.isclose(segments[0].pPoint_start.dLatitude_degree, 0.0)
        assert np.isclose(segments[-1].pPoint_end.dLatitude_degree, 10.0)


class TestPyCircle:
    def test_circle_init(self):
        center = pypoint({"dLongitude_degree": 10.0, "dLatitude_degree": 20.0})
        circle_pts = [
            pypoint({"dLongitude_degree": 10.1, "dLatitude_degree": 20.0}),
            pypoint({"dLongitude_degree": 9.9, "dLatitude_degree": 20.0}),
        ]
        circle = pycircle(center, circle_pts)
        assert circle.pPoint_center == center
        assert len(circle.aPoint_circle) == 2
        js = circle.tojson()
        assert "pPoint_center" in js
