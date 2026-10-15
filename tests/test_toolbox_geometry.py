import numpy as np
import pytest

from pyearth.toolbox.geometry import calculate_hexagon_area


def test_calculate_hexagon_area_scalar():
    area = calculate_hexagon_area(1.0)
    expected = 1.5 * np.sqrt(3.0)
    assert np.isclose(area, expected)

    assert calculate_hexagon_area(0.0) == 0.0


def test_calculate_hexagon_area_array():
    edges = np.array([1.0, 2.0, 3.0])
    areas = calculate_hexagon_area(edges)
    expected = 1.5 * np.sqrt(3.0) * edges**2
    assert np.allclose(areas, expected)


def test_calculate_hexagon_area_invalid():
    with pytest.raises(TypeError):
        calculate_hexagon_area(None)
    with pytest.raises(ValueError):
        calculate_hexagon_area(-1.0)
