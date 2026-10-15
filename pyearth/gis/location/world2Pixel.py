"""Convert between geographic/projected coordinates and raster pixel coordinates."""

from typing import List, Sequence, Tuple, Union

__all__ = [
    "world_to_pixel",
    "pixel_to_world",
]


def world_to_pixel(
    geotransform: Sequence[float], x: float, y: float
) -> Tuple[int, int]:
    """Convert world coordinates to raster pixel indices using a GDAL geotransform.

    Parameters
    ----------
    geotransform : Sequence[float]
        GDAL geotransform tuple/list with 6 elements:
        - [0]: X origin (top-left longitude or easting)
        - [1]: Pixel width (W-E resolution)
        - [2]: Row rotation (typically 0 for north-up images)
        - [3]: Y origin (top-left latitude or northing)
        - [4]: Column rotation (typically 0 for north-up images)
        - [5]: Pixel height (N-S resolution, typically negative)
    x : float
        X coordinate (longitude or easting) in the raster's CRS.
    y : float
        Y coordinate (latitude or northing) in the raster's CRS.

    Returns
    -------
    Tuple[int, int]
        (column_index, row_index) 0-based pixel indices from top-left.

    Raises
    ------
    ValueError
        If geotransform doesn't have 6 elements or determinant is zero.
    TypeError
        If coordinates or geotransform elements are non-numeric.
    """
    if len(geotransform) != 6:
        raise ValueError(
            f"Geotransform must have exactly 6 elements. Got {len(geotransform)}."
        )

    try:
        x_val = float(x)
        y_val = float(y)
        c0, c1, c2, c3, c4, c5 = [float(v) for v in geotransform]
    except (TypeError, ValueError) as e:
        raise TypeError(f"Coordinates and geotransform values must be numeric: {e}")

    # Standard north-up fast path
    if c2 == 0.0 and c4 == 0.0:
        if c1 == 0.0 or c5 == 0.0:
            raise ValueError(
                "Pixel dimensions cannot be zero. Invalid geotransform."
            )
        col = int((x_val - c0) / c1)
        row = int((y_val - c3) / c5)
        return col, row

    # General affine transformation with rotation
    det = c1 * c5 - c2 * c4
    if det == 0.0:
        raise ValueError(
            "Geotransform determinant is zero; transformation is not invertible."
        )

    dx = x_val - c0
    dy = y_val - c3
    col = int((dx * c5 - dy * c2) / det)
    row = int((dy * c1 - dx * c4) / det)

    return col, row


def pixel_to_world(
    geotransform: Sequence[float],
    col: float,
    row: float,
    center: bool = False,
) -> Tuple[float, float]:
    """Convert raster pixel indices to world coordinates using a GDAL geotransform.

    Parameters
    ----------
    geotransform : Sequence[float]
        GDAL geotransform tuple/list with 6 elements.
    col : float
        Pixel column index (0-based from left).
    row : float
        Pixel row index (0-based from top).
    center : bool, optional
        If True, return coordinates of the pixel center (offset by +0.5).
        If False (default), return coordinates of the pixel top-left corner.

    Returns
    -------
    Tuple[float, float]
        (x, y) world coordinates.
    """
    if len(geotransform) != 6:
        raise ValueError(
            f"Geotransform must have exactly 6 elements. Got {len(geotransform)}."
        )

    try:
        c0, c1, c2, c3, c4, c5 = [float(v) for v in geotransform]
        c = float(col) + (0.5 if center else 0.0)
        r = float(row) + (0.5 if center else 0.0)
    except (TypeError, ValueError) as e:
        raise TypeError(f"Pixel indices and geotransform must be numeric: {e}")

    x = c0 + c * c1 + r * c2
    y = c3 + c * c4 + r * c5
    return x, y


