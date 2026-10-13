import numpy as np
from pyearth.gis.idl_handler import (
    split_international_date_line_polygon_coordinates,
)


def test_idl_split():
    # Example 1: Simple rectangular polygon crossing IDL
    coords = np.array([
        [170.0, -10.0],
        [170.0,  10.0],
        [-170.0, 10.0],
        [-170.0, -10.0],
        [170.0, -10.0],
    ])
    eastern, western = split_international_date_line_polygon_coordinates(coords)
    assert np.all(eastern[:, 0] > 0)
    assert np.all(western[:, 0] < 0)
    assert np.allclose(eastern[0], eastern[-1])
    assert np.allclose(western[0], western[-1])

    # Example 4: Complex polygon
    coords_complex = np.array([
        [175.0, 0.0],
        [178.0, 2.0],
        [179.0, 5.0],
        [-179.0, 5.0],
        [-178.0, 2.0],
        [-175.0, 0.0],
        [-178.0, -2.0],
        [-179.0, -5.0],
        [179.0, -5.0],
        [178.0, -2.0],
        [175.0, 0.0],
    ])
    eastern2, western2 = split_international_date_line_polygon_coordinates(coords_complex)
    assert np.all(eastern2[:, 0] > 0)
    assert np.all(western2[:, 0] < 0)
    assert np.allclose(eastern2[0], eastern2[-1])
    assert np.allclose(western2[0], western2[-1])

    # IDL-vertex-only path: polygon touching ±180 but not crossing
    coords_touch = np.array([
        [170.0, 0.0],
        [180.0, 0.0],
        [180.0, 10.0],
        [170.0, 10.0],
        [170.0, 0.0],
    ])
    result_touch = split_international_date_line_polygon_coordinates(coords_touch)
    assert result_touch[0].shape[0] > 0
    assert result_touch[1].shape[0] == 0


if __name__ == "__main__":
    test_idl_split()
    print("All IDL split tests passed successfully.")
