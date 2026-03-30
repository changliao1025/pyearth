import numpy as np
from typing import Union
from pyearth.toolbox.data.remove_duplicate_closure import remove_duplicate_closure as remove_duplicate_closure

def calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    x, y = coords[:, 0], coords[:, 1]

    # Vectorized shoelace formula - more efficient than loops
    # Handle the wrap-around (last point to first point) implicitly
    x_rolled = np.roll(x, -1)  # x[i+1] for all i, with wrap-around
    y_rolled = np.roll(y, -1)  # y[i+1] for all i, with wrap-around

    signed_area = 0.5 * np.sum(x * y_rolled - x_rolled * y)
    return signed_area


def calculate_signed_area_spherical_polar(
    coords: np.ndarray, pole: str = "north"
) -> float:
    """Calculate signed spherical area for polygons enclosing a pole.

    The polygon is projected with Lambert azimuthal equal-area (LAEA) centered
    on the requested pole, then planar signed area is computed with the
    shoelace formula. With unit sphere radius, output area is in steradians.

    Parameters
    ----------
    coords : np.ndarray
        Polygon coordinates as (lon, lat) degrees.
    pole : str, default="north"
        Polar center for projection, either "north" or "south".

    Returns
    -------
    float
        Signed spherical area (steradians on unit sphere).
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return 0.0

    arr = remove_duplicate_closure(coords)
    if len(arr) < 3:
        return 0.0

    lon_rad = np.deg2rad(arr[:, 0])
    lat_rad = np.deg2rad(arr[:, 1])

    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)

    pole_lc = pole.lower()
    if pole_lc == "north":
        # LAEA centered at +90 deg latitude
        denom = 1.0 + sin_lat
        # Avoid divide-by-zero at antipode (not expected for pole-enclosing cells)
        denom = np.maximum(denom, 1.0e-15)
        k = np.sqrt(2.0 / denom)
        x = k * cos_lat * np.sin(lon_rad)
        y = -k * cos_lat * np.cos(lon_rad)
    elif pole_lc == "south":
        # LAEA centered at -90 deg latitude
        denom = 1.0 - sin_lat
        denom = np.maximum(denom, 1.0e-15)
        k = np.sqrt(2.0 / denom)
        x = k * cos_lat * np.sin(lon_rad)
        y = k * cos_lat * np.cos(lon_rad)
    else:
        raise ValueError("pole must be either 'north' or 'south'")

    projected = np.column_stack((x, y))
    return calculate_signed_area_shoelace(projected)


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
    from pyearth.gis.geometry.pole_check import polygon_includes_pole

    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    if len(coords) < 3:
        return True  # Degenerate case

    if polygon_includes_pole(coords, pole="north"):
        signed_area = calculate_signed_area_spherical_polar(coords, pole="north")
    elif polygon_includes_pole(coords, pole="south"):
        signed_area = calculate_signed_area_spherical_polar(coords, pole="south")
    else:
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
