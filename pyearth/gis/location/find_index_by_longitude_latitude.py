"""Grid indexing utilities for finding cell indices from geographic coordinates."""

from typing import Tuple

__all__ = [
    "find_index_by_longitude_latitude",
]


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
    2D grid/raster. Handles both cell-center and corner registration modes.

    Parameters
    ----------
    dLon_in : float
        Longitude of the target point in decimal degrees.
    dLat_in : float
        Latitude of the target point in decimal degrees.
    dLon_min_in : float
        Minimum longitude (western edge) of the grid in decimal degrees.
    dLat_max_in : float
        Maximum latitude (northern edge) of the grid in decimal degrees.
    dResolution_x_in : float
        Grid cell width in the longitude (X) direction, in degrees. Must be positive.
    dResolution_y_in : float
        Grid cell height in the latitude (Y) direction, in degrees. Must be positive.
    iFlag_center_in : bool, optional
        Grid registration mode. Default is True.
        - True: Cell-center registration (dLon_min/dLat_max are cell centers).
        - False: Corner registration (dLon_min/dLat_max are upper-left corners).

    Returns
    -------
    Tuple[int, int]
        Grid indices (row_index, column_index) where:
        - row_index: 0-based index counting from north to south.
        - column_index: 0-based index counting from west to east.

    Raises
    ------
    ValueError
        If resolutions are not positive or coordinates are out of bounds.
    TypeError
        If input values are non-numeric.
    """
    try:
        lon = float(dLon_in)
        lat = float(dLat_in)
        lon_min = float(dLon_min_in)
        lat_max = float(dLat_max_in)
        res_x = float(dResolution_x_in)
        res_y = float(dResolution_y_in)
    except (TypeError, ValueError) as e:
        raise TypeError(f"All coordinate and resolution values must be numeric: {e}")

    if res_x <= 0:
        raise ValueError(f"X resolution must be positive. Got {res_x}")
    if res_y <= 0:
        raise ValueError(f"Y resolution must be positive. Got {res_y}")

    if iFlag_center_in:
        row_index = round((lat_max - lat) / res_y)
        column_index = round((lon - lon_min) / res_x)
    else:
        lat_max_center = lat_max - 0.5 * res_y
        lon_min_center = lon_min + 0.5 * res_x
        row_index = round((lat_max_center - lat) / res_y)
        column_index = round((lon - lon_min_center) / res_x)

    if row_index < 0:
        raise ValueError(
            f"Point latitude {lat}° is north of grid bounds (lat_max={lat_max}°)."
        )
    if column_index < 0:
        raise ValueError(
            f"Point longitude {lon}° is west of grid bounds (lon_min={lon_min}°)."
        )

    return int(row_index), int(column_index)
