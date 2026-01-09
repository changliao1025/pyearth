"""
Grid indexing utilities for finding cell indices from geographic coordinates.

This module provides functions for converting longitude/latitude coordinates
to grid cell indices (row, column) in regular 2D raster grids.
"""

from typing import Tuple


def find_index_by_longitude_latitude(
    dLon_in: float,
    dLat_in: float,
    dLon_min_in: float,
    dLat_max_in: float,
    dResolution_x_in: float,
    dResolution_y_in: float,
    iFlag_center_in: bool = True,
) -> Tuple[int, int]:
    """Find the grid cell indices (row, column) for a geographic coordinate.

    Converts a longitude/latitude point to grid cell indices within a regular
    2D grid/raster. The function handles both cell-center and corner registration
    systems commonly used in raster datasets.

    Parameters
    ----------
    dLon_in : float
        Longitude of the target point in decimal degrees.
        Should be in the same range as the grid (either [-180, 180] or [0, 360]).
    dLat_in : float
        Latitude of the target point in decimal degrees.
        Valid range: [-90, 90].
    dLon_min_in : float
        Minimum longitude (western edge) of the grid in decimal degrees.
        For cell-center registration, this is the longitude of the westernmost cell center.
        For corner registration, this is the longitude of the western grid boundary.
    dLat_max_in : float
        Maximum latitude (northern edge) of the grid in decimal degrees.
        For cell-center registration, this is the latitude of the northernmost cell center.
        For corner registration, this is the latitude of the northern grid boundary.
    dResolution_x_in : float
        Grid cell width in the longitude (X) direction, in degrees.
        Must be positive.
    dResolution_y_in : float
        Grid cell height in the latitude (Y) direction, in degrees.
        Must be positive.
    iFlag_center_in : bool, optional
        Grid registration mode. Default is True.
        - True: Cell-center registration. dLon_min/dLat_max are cell centers.
        - False: Corner registration. dLon_min/dLat_max are upper-left corners.

    Returns
    -------
    Tuple[int, int]
        Grid indices (row_index, column_index) where:
        - row_index: Row index (0-based, counting from north to south)
        - column_index: Column index (0-based, counting from west to east)

    Raises
    ------
    ValueError
        If resolutions are not positive.
        If the point is outside the expected grid bounds.
    TypeError
        If input values cannot be converted to float.

    Notes
    -----
    - Row indices increase from north to south (top to bottom)
    - Column indices increase from west to east (left to right)
    - The function uses rounding, so points near cell boundaries may round
      to either adjacent cell
    - For corner registration, the function adjusts to find the cell center
      before calculating indices

    Grid Registration Modes
    ------------------------
    Cell-center registration:
        Grid coordinates represent the centers of grid cells.
        dLon_min is at the center of the leftmost column.
        dLat_max is at the center of the topmost row.

    Corner registration:
        Grid coordinates represent the upper-left corners of grid cells.
        dLon_min is at the left edge of the leftmost column.
        dLat_max is at the top edge of the topmost row.

    Examples
    --------
    >>> # Cell-center registration: 1° grid starting at (0°E, 50°N)
    >>> row, col = find_index_by_longitude_latitude(
    ...     dLon_in=3.5, dLat_in=47.5,
    ...     dLon_min_in=0.0, dLat_max_in=50.0,
    ...     dResolution_x_in=1.0, dResolution_y_in=1.0,
    ...     iFlag_center_in=True
    ... )
    >>> row, col
    (2, 4)

    >>> # Corner registration: same grid but corners defined
    >>> row, col = find_index_by_longitude_latitude(
    ...     dLon_in=3.5, dLat_in=47.5,
    ...     dLon_min_in=-0.5, dLat_max_in=50.5,
    ...     dResolution_x_in=1.0, dResolution_y_in=1.0,
    ...     iFlag_center_in=False
    ... )
    >>> row, col
    (2, 4)

    See Also
    --------
    numpy.digitize : For more complex binning operations
    """
    # Validate and convert inputs
    try:
        lon = float(dLon_in)
        lat = float(dLat_in)
        lon_min = float(dLon_min_in)
        lat_max = float(dLat_max_in)
        res_x = float(dResolution_x_in)
        res_y = float(dResolution_y_in)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"All coordinate and resolution values must be numeric. Error: {e}"
        )

    # Validate resolutions
    if res_x <= 0:
        raise ValueError(f"X resolution must be positive. Got {res_x}")
    if res_y <= 0:
        raise ValueError(f"Y resolution must be positive. Got {res_y}")

    # Calculate indices based on registration mode
    if iFlag_center_in:
        # Cell-center registration: min/max are at cell centers
        # Row index: how many cells down from the top (north to south)
        row_index = round((lat_max - lat) / res_y)

        # Column index: how many cells right from the left (west to east)
        column_index = round((lon - lon_min) / res_x)
    else:
        # Corner registration: min/max are at upper-left corners
        # Adjust to cell centers before calculating
        # Shift by half a cell to convert corner to center reference
        lat_max_center = lat_max - 0.5 * res_y
        lon_min_center = lon_min + 0.5 * res_x

        row_index = round((lat_max_center - lat) / res_y)
        column_index = round((lon - lon_min_center) / res_x)

    # Validate that indices are non-negative (basic bounds check)
    if row_index < 0:
        raise ValueError(
            f"Point latitude {lat}° is north of grid bounds (lat_max={lat_max}°). "
            f"Calculated row index: {row_index}"
        )
    if column_index < 0:
        raise ValueError(
            f"Point longitude {lon}° is west of grid bounds (lon_min={lon_min}°). "
            f"Calculated column index: {column_index}"
        )

    return int(row_index), int(column_index)
