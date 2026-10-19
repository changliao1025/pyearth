import numpy as np

from pyearth.gis.idl_handler import IdlHandler, IdlStrategy


def test_idl_split():
    handler = IdlHandler(strategy=IdlStrategy.SPLIT)

    # Example 1: Simple rectangular polygon crossing IDL
    coords = np.array(
        [
            [170.0, -10.0],
            [170.0, 10.0],
            [-170.0, 10.0],
            [-170.0, -10.0],
            [170.0, -10.0],
        ]
    )
    result = handler.detect(coords)
    assert result.crosses_idl is True
    assert result.sub_polygons is not None
    eastern, western = result.sub_polygons
    assert np.all(eastern[:, 0] > 0)
    assert np.all(western[:, 0] < 0)
    assert np.allclose(eastern[0], eastern[-1])
    assert np.allclose(western[0], western[-1])

    # Example 4: Complex polygon
    coords_complex = np.array(
        [
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
        ]
    )
    result2 = handler.detect(coords_complex)
    assert result2.crosses_idl is True
    assert result2.sub_polygons is not None
    eastern2, western2 = result2.sub_polygons
    assert np.all(eastern2[:, 0] > 0)
    assert np.all(western2[:, 0] < 0)
    assert np.allclose(eastern2[0], eastern2[-1])
    assert np.allclose(western2[0], western2[-1])

    # IDL-vertex-only path: polygon touching ±180 but not crossing
    coords_touch = np.array(
        [
            [170.0, 0.0],
            [180.0, 0.0],
            [180.0, 10.0],
            [170.0, 10.0],
            [170.0, 0.0],
        ]
    )
    result_touch = handler.detect(coords_touch)
    assert result_touch.crosses_idl is False
    assert result_touch.adjusted_coords is not None


if __name__ == "__main__":
    test_idl_split()
    print("All IDL split tests passed successfully.")
