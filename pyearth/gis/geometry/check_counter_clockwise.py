import numpy as np
from typing import Union


def calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    x, y = coords[:, 0], coords[:, 1]

    # Vectorized shoelace formula - more efficient than loops
    # Handle the wrap-around (last point to first point) implicitly
    x_rolled = np.roll(x, -1)  # x[i+1] for all i, with wrap-around
    y_rolled = np.roll(y, -1)  # y[i+1] for all i, with wrap-around

    signed_area = 0.5 * np.sum(x * y_rolled - x_rolled * y)
    return signed_area


def check_counter_clockwise(coords: np.ndarray) -> bool:
    """Check if polygon coordinates are in counter-clockwise order.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    bool
        True if vertices are in counter-clockwise order, False otherwise.
    """
    from pyearth.gis.geometry.international_date_line_utility import (
        check_cross_international_date_line_polygon,
        unwrap_longitudes,
    )

    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return True  # Degenerate case

    # Check if polygon crosses the International Date Line
    iFlag_cross, _ = check_cross_international_date_line_polygon(coords)
    if iFlag_cross:
        coords_unwrapped = unwrap_longitudes(coords)
        # Calculate signed area using optimized shoelace formula
        signed_area = calculate_signed_area_shoelace(coords_unwrapped)
    else:
        # Standard case: calculate signed area directly
        signed_area = calculate_signed_area_shoelace(coords)
    return signed_area > 0


def get_polygon_orientation(coords: np.ndarray) -> str:
    """Get the orientation of polygon coordinates as a string.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    str
        Either "counter-clockwise" or "clockwise".

    Examples
    --------
    >>> coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
    >>> get_polygon_orientation(coords)
    'counter-clockwise'
    """
    return "counter-clockwise" if check_counter_clockwise(coords) else "clockwise"
