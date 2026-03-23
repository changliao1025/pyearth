import numpy as np
from pyearth.gis.geometry.international_date_line_utility import (
    split_international_date_line_polygon_coordinates,
)

# Example 1: Simple rectangular polygon crossing IDL
coords = np.array([
    [170.0, -10.0],
    [170.0,  10.0],
    [-170.0, 10.0],
    [-170.0,-10.0],
    [170.0, -10.0],
])
eastern, western = split_international_date_line_polygon_coordinates(coords)
print('=== Example 1 ===')
print('eastern:\n', eastern)
print('western:\n', western)
print('all eastern lon > 0:', np.all(eastern[:, 0] > 0))
print('all western lon < 0:', np.all(western[:, 0] < 0))
print('eastern closed:', np.allclose(eastern[0], eastern[-1]))
print('western closed:', np.allclose(western[0], western[-1]))

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
    [175.0, 0.0]
])
eastern2, western2 = split_international_date_line_polygon_coordinates(coords_complex)
print()
print('=== Example 4 ===')
print('eastern:\n', eastern2)
print('western:\n', western2)
print('all eastern lon > 0:', np.all(eastern2[:, 0] > 0))
print('all western lon < 0:', np.all(western2[:, 0] < 0))
print('eastern closed:', np.allclose(eastern2[0], eastern2[-1]))
print('western closed:', np.allclose(western2[0], western2[-1]))

# IDL-vertex-only path: polygon touching ±180 but not crossing
coords_touch = np.array([
    [170.0, 0.0],
    [180.0, 0.0],
    [180.0, 10.0],
    [170.0, 10.0],
    [170.0, 0.0],
])
result_touch = split_international_date_line_polygon_coordinates(coords_touch)
print()
print('=== IDL-vertex-only (touch, no cross) ===')
print('result[0] shape:', result_touch[0].shape)
print('result[1] shape:', result_touch[1].shape)
print('result[0]:\n', result_touch[0])
