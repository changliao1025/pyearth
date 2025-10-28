import numpy as np
from typing import Union


def unwrap_longitudes(coords: np.ndarray) -> np.ndarray:
    """Unwrap longitude coordinates to handle International Date Line crossings.

    Adjusts longitude values to ensure they are all within 180 degrees of the
    first coordinate, preventing artificial jumps when calculating polygon areas.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    np.ndarray
        Array with unwrapped longitude coordinates, same shape as input.

    Notes
    -----
    This function uses the first longitude point as a reference and adjusts all
    subsequent longitudes to be within ±180 degrees of this reference point by
    adding or subtracting 360 degrees as needed.

    Examples
    --------
    >>> coords = np.array([[-170, 10], [170, 20], [-160, 30]])
    >>> unwrapped = unwrap_longitudes(coords)
    >>> # Longitudes adjusted to avoid large jumps across IDL
    """
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")

    coords_unwrapped = coords.copy()
    lons = coords_unwrapped[:, 0]
    ref_lon = lons[0]

    # Vectorized approach for better performance
    diff = lons - ref_lon
    # Adjust longitudes that are more than 180 degrees away
    lons[diff > 180] -= 360
    lons[diff < -180] += 360

    return coords_unwrapped


def _crosses_international_date_line(coords: np.ndarray) -> bool:
    """Check if polygon coordinates cross the International Date Line.

    Parameters
    ----------
    coords : np.ndarray
        Array of polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    bool
        True if the polygon crosses the IDL, False otherwise.

    Notes
    -----
    Uses a more robust method than simple longitude range checking by
    examining consecutive longitude differences.
    """
    lons = coords[:, 0]

    # Check for large jumps between consecutive points (> 180 degrees)
    lon_diffs = np.abs(np.diff(lons))
    max_jump = np.max(lon_diffs)

    # Also check the wrap-around from last to first point
    wrap_jump = abs(lons[-1] - lons[0])

    return max_jump > 180 or wrap_jump > 180


def _calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    """Calculate signed area using the shoelace formula (optimized version).

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    float
        Signed area of the polygon. Positive for CCW, negative for CW.

    Notes
    -----
    Uses vectorized operations for optimal performance.
    Formula: Area = 0.5 * sum(x[i]*y[i+1] - x[i+1]*y[i])
    """
    x, y = coords[:, 0], coords[:, 1]

    # Vectorized shoelace formula - more efficient than loops
    # Handle the wrap-around (last point to first point) implicitly
    x_rolled = np.roll(x, -1)  # x[i+1] for all i, with wrap-around
    y_rolled = np.roll(y, -1)  # y[i+1] for all i, with wrap-around

    signed_area = 0.5 * np.sum(x * y_rolled - x_rolled * y)
    return signed_area


def check_ccw_idl(coords: np.ndarray) -> bool:
    """Check counter-clockwise orientation for polygons crossing the International Date Line.

    This function specifically handles polygons that cross the International Date Line
    by unwrapping longitude coordinates before applying the standard CCW check.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    bool
        True if coordinates are in counter-clockwise order, False if clockwise.

    Raises
    ------
    ValueError
        If coords is not a valid polygon coordinate array.

    Notes
    -----
    This function is automatically called by check_ccw() when IDL crossing is detected.
    It unwraps longitudes to handle the coordinate discontinuity at ±180 degrees.

    Examples
    --------
    >>> # Polygon crossing IDL in CCW order
    >>> coords = np.array([[170, 10], [-170, 20], [-160, 30], [160, 40], [170, 10]])
    >>> check_ccw_idl(coords)
    True
    """
    # Validate input using the same validation as check_ccw
    _validate_polygon_coords(coords)

    # Unwrap longitudes to handle IDL crossing
    coords_unwrapped = unwrap_longitudes(coords)

    # Calculate signed area using optimized shoelace formula
    signed_area = _calculate_signed_area_shoelace(coords_unwrapped)

    return signed_area > 0


def _validate_polygon_coords(coords: np.ndarray) -> None:
    """Validate polygon coordinate array.

    Parameters
    ----------
    coords : np.ndarray
        Array to validate.

    Raises
    ------
    ValueError
        If coords is not a valid polygon coordinate array.
    """
    if not isinstance(coords, np.ndarray):
        raise ValueError("Coordinates must be a numpy array")

    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(
            f"Coordinates must be a 2D array with shape (n, 2). Got shape {coords.shape}"
        )

    if len(coords) < 3:
        raise ValueError(
            f"Polygon must have at least 3 points. Got {len(coords)} points"
        )


def check_ccw(coords: np.ndarray) -> bool:
    """Determine if polygon coordinates are in counter-clockwise (CCW) order.

    This function automatically detects and handles polygons that cross the
    International Date Line (IDL) by using appropriate coordinate unwrapping.
    Uses the optimized shoelace formula for area calculation.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates in (longitude, latitude) format.

    Returns
    -------
    bool
        True if coordinates are in counter-clockwise order, False if clockwise.

    Raises
    ------
    ValueError
        If coords is not a 2D array with shape (n, 2) where n >= 3.

    Notes
    -----
    - Uses the shoelace formula: Area = 0.5 * sum(x[i]*y[i+1] - x[i+1]*y[i])
    - Positive signed area indicates CCW, negative indicates CW
    - Automatically handles International Date Line crossings
    - Optimized with vectorized operations for better performance
    - Minimum 3 points required for a valid polygon

    Examples
    --------
    >>> # Square in CCW order
    >>> coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
    >>> check_ccw(coords)
    True

    >>> # Square in CW order
    >>> coords = np.array([[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]])
    >>> check_ccw(coords)
    False

    >>> # Polygon crossing IDL in CCW order
    >>> coords = np.array([[170, 10], [-170, 20], [-160, 30], [160, 40], [170, 10]])
    >>> check_ccw(coords)
    True
    """
    # Validate input coordinates
    _validate_polygon_coords(coords)

    # Check if polygon crosses the International Date Line
    if _crosses_international_date_line(coords):
        # Use specialized IDL handling
        return check_ccw_idl(coords)
    else:
        # Standard case: calculate signed area directly
        signed_area = _calculate_signed_area_shoelace(coords)
        return signed_area > 0


def is_polygon_clockwise(coords: np.ndarray) -> bool:
    """Determine if polygon coordinates are in clockwise (CW) order.

    Convenience function that returns the opposite of check_ccw().

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (n, 2) representing polygon coordinates.

    Returns
    -------
    bool
        True if coordinates are in clockwise order, False if counter-clockwise.

    Examples
    --------
    >>> coords = np.array([[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]])
    >>> is_polygon_clockwise(coords)
    True
    """
    return not check_ccw(coords)


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
    return "counter-clockwise" if check_ccw(coords) else "clockwise"