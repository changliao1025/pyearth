import numpy as np

from pyearth.gis.polar_handler import PolarHandler, polygon_includes_pole


def test_polygon_includes_pole_for_south_pole_enclosing_polygon():
    coords = np.array(
        [
            [0.0, -80.0],
            [90.0, -80.0],
            [180.0, -80.0],
            [-90.0, -80.0],
            [0.0, -80.0],
        ]
    )

    assert polygon_includes_pole(coords, pole="south") is True
    assert polygon_includes_pole(coords, pole="north") is False


def test_polar_handler_detect_and_includes_match_wrapper():
    coords = np.array(
        [
            [0.0, -80.0],
            [90.0, -80.0],
            [180.0, -80.0],
            [-90.0, -80.0],
            [0.0, -80.0],
        ]
    )

    handler = PolarHandler(pole="south")
    result = handler.detect(coords)

    assert result.includes_pole is True
    assert handler.includes(coords) is True
    assert result.projected_coords is not None
    assert result.projected_coords.shape[1] == 2
