"""
Convert geographic/projected coordinates to raster pixel coordinates.

This module provides utilities for converting world coordinates (geographic or projected)
to pixel/line indices in a georeferenced raster using GDAL geotransform parameters.
"""

from typing import Tuple, Union, List


def world2Pixel(
    pGeoMatrix_in: Union[Tuple, List], dx_in: float, dy_in: float
) -> Tuple[int, int]:
    """Convert world coordinates to pixel indices using a GDAL geotransform.

    Transforms geographic or projected coordinates to raster pixel/line indices
    using the affine transformation parameters from a GDAL geotransform matrix.

    Parameters
    ----------
    pGeoMatrix_in : tuple or list
        GDAL geotransform tuple with 6 elements, typically from
        dataset.GetGeoTransform(). Format: (x_origin, pixel_width, x_rotation,
        y_origin, y_rotation, pixel_height).
        - [0]: X coordinate of upper-left corner (origin)
        - [1]: Pixel width (W-E pixel resolution)
        - [2]: Row rotation (0 for north-up images)
        - [3]: Y coordinate of upper-left corner (origin)
        - [4]: Column rotation (0 for north-up images)
        - [5]: Pixel height (N-S pixel resolution, typically negative)
    dx_in : float
        X coordinate (longitude or easting) in the raster's spatial reference system.
    dy_in : float
        Y coordinate (latitude or northing) in the raster's spatial reference system.

    Returns
    -------
    Tuple[int, int]
        Pixel indices (column, row) where:
        - column: Pixel index in X direction (0-based, from left)
        - row: Pixel index in Y direction (0-based, from top)

    Raises
    ------
    ValueError
        If geotransform doesn't have exactly 6 elements.
        If pixel dimensions are zero (invalid geotransform).
    TypeError
        If inputs cannot be converted to appropriate types.

    Notes
    -----
    - For north-up images without rotation (most common):
      * pixel_col = (x - x_origin) / pixel_width
      * pixel_row = (y_origin - y) / |pixel_height|
    - Pixel height ([5]) is typically negative for north-up images
    - Rotation parameters ([2], [4]) are usually 0 for north-up images
    - This function assumes no rotation; rotation parameters are extracted but not used
    - Coordinates outside raster bounds will return negative or large indices
    - Results are truncated to integers (floor division)

    GDAL Geotransform Format
    -------------------------
    The geotransform defines the relationship between raster pixel/line coordinates
    and georeferenced coordinates:
        X_geo = geotransform[0] + pixel*geotransform[1] + line*geotransform[2]
        Y_geo = geotransform[3] + pixel*geotransform[4] + line*geotransform[5]

    This function inverts this transformation to find pixel/line from X_geo/Y_geo.

    Examples
    --------
    >>> # Geotransform for 1-degree grid starting at (0°E, 90°N)
    >>> geotransform = (0.0, 1.0, 0.0, 90.0, 0.0, -1.0)
    >>> # Point at (45°E, 45°N)
    >>> col, row = world2Pixel(geotransform, 45.0, 45.0)
    >>> col, row
    (45, 45)

    >>> # Geotransform for UTM zone with 30m pixels
    >>> geotransform = (500000.0, 30.0, 0.0, 4500000.0, 0.0, -30.0)
    >>> # Point 300m east and 600m south of origin
    >>> col, row = world2Pixel(geotransform, 500300.0, 4499400.0)
    >>> col, row
    (10, 20)

    See Also
    --------
    gdal.Dataset.GetGeoTransform : Get geotransform from GDAL dataset

    References
    ----------
    - GDAL Geotransform Tutorial:
      https://gdal.org/tutorials/geotransforms_tut.html
    """
    # Validate geotransform
    if not isinstance(pGeoMatrix_in, (tuple, list)):
        raise TypeError(
            f"Geotransform must be a tuple or list. Got {type(pGeoMatrix_in).__name__}"
        )

    if len(pGeoMatrix_in) != 6:
        raise ValueError(
            f"Geotransform must have exactly 6 elements. Got {len(pGeoMatrix_in)}"
        )

    # Validate and convert coordinates
    try:
        x_coord = float(dx_in)
        y_coord = float(dy_in)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Coordinates must be numeric values. "
            f"Got x={dx_in}, y={dy_in}. Error: {e}"
        )

    # Extract geotransform parameters
    x_origin = pGeoMatrix_in[0]  # Upper-left X (longitude/easting)
    pixel_width = pGeoMatrix_in[1]  # W-E pixel resolution (positive)
    x_rotation = pGeoMatrix_in[2]  # Row rotation (typically 0)
    y_origin = pGeoMatrix_in[3]  # Upper-left Y (latitude/northing)
    y_rotation = pGeoMatrix_in[4]  # Column rotation (typically 0)
    pixel_height = pGeoMatrix_in[5]  # N-S pixel resolution (typically negative)

    # Validate pixel dimensions are non-zero
    if pixel_width == 0:
        raise ValueError(
            "Pixel width (geotransform[1]) cannot be zero. Invalid geotransform."
        )
    if pixel_height == 0:
        raise ValueError(
            "Pixel height (geotransform[5]) cannot be zero. Invalid geotransform."
        )

    # Calculate pixel column (X direction)
    # For north-up images: pixel = (x - x_origin) / pixel_width
    pixel_col = int((x_coord - x_origin) / pixel_width)

    # Calculate pixel row (Y direction)
    # For north-up images with negative pixel_height:
    # row = (y_origin - y) / abs(pixel_height)
    # This correctly uses pixel_height (not pixel_width as in the original bug)
    pixel_row = int((y_origin - y_coord) / abs(pixel_height))

    return pixel_col, pixel_row
