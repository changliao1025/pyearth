"""Unified polar region handler for polygon processing.

This module provides a centralized approach to detecting whether polygons
enclose the North or South Pole.  It projects lon/lat vertices into a local
polar plane where the pole maps to the origin, then runs a standard 2D
point-in-polygon test.

Usage:
    from pyearth.gis.polar_handler import PolarHandler, PoleType

    handler = PolarHandler(pole=PoleType.SOUTH)

    # Detect whether a polygon includes the pole
    result = handler.detect(coords_array)
    print(result.includes_pole)

    # Quick functional check
    if handler.includes(coords_array):
        ...

    # Or use the module-level convenience function
    from pyearth.gis.polar_handler import polygon_includes_pole
    polygon_includes_pole(coords_array, pole="south")
"""


from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from pyearth.toolbox.data.remove_duplicate_closure import (
    remove_duplicate_closure as remove_duplicate_closure,
)


# ---------------------------------------------------------------------------
# Module-level constants (single source of truth)
# ---------------------------------------------------------------------------
POLAR_TOLERANCE = 1.0e-12  # Tolerance for point-on-segment / point-in-polygon tests
LATITUDE_BOUNDARY = 90.0   # Absolute latitude of the poles


# ---------------------------------------------------------------------------
# Enums & Data classes
# ---------------------------------------------------------------------------
class PoleType(Enum):
    """Which pole to test against.

    Attributes
    ----------
    NORTH : str
        Test whether the polygon encloses the North Pole (90°N).
    SOUTH : str
        Test whether the polygon encloses the South Pole (90°S).
    """

    NORTH = "north"
    SOUTH = "south"


@dataclass
class PolarResult:
    """Encapsulates the result of a polar-inclusion test.

    Attributes
    ----------
    includes_pole : bool
        Whether the input polygon encloses the requested pole.
    pole : PoleType
        The pole that was tested.
    projected_coords : numpy.ndarray or None
        The 2-D projected coordinates used for the point-in-polygon test
        (the pole maps to the origin).  ``None`` when the input was invalid.
    """

    includes_pole: bool
    pole: PoleType
    projected_coords: Optional[np.ndarray] = None


# ---------------------------------------------------------------------------
# Helpers: 2-D geometry primitives
# ---------------------------------------------------------------------------
def _point_on_segment_2d(
    point: NDArray[np.floating],
    seg_start: NDArray[np.floating],
    seg_end: NDArray[np.floating],
    tol: float = POLAR_TOLERANCE,
) -> bool:
    """Return True if a 2-D point lies on a line segment within tolerance.

    Parameters
    ----------
    point : array_like, shape (2,)
        The query point.
    seg_start, seg_end : array_like, shape (2,)
        Endpoints of the line segment.
    tol : float
        Numerical tolerance.

    Returns
    -------
    bool
    """
    px, py = point
    x1, y1 = seg_start
    x2, y2 = seg_end

    dx = x2 - x1
    dy = y2 - y1

    # Degenerate segment
    if abs(dx) < tol and abs(dy) < tol:
        return np.hypot(px - x1, py - y1) <= tol

    # Cross-product distance to the infinite line
    cross = (px - x1) * dy - (py - y1) * dx
    if abs(cross) > tol:
        return False

    # Dot-product bounds check for segment extents
    dot = (px - x1) * dx + (py - y1) * dy
    if dot < -tol:
        return False

    seg_len_sq = dx * dx + dy * dy
    if dot - seg_len_sq > tol:
        return False

    return True


def _point_in_polygon_2d(
    point: NDArray[np.floating],
    polygon: NDArray[np.floating],
    include_boundary: bool = False,
    tol: float = POLAR_TOLERANCE,
) -> bool:
    """2-D point-in-polygon using ray casting.

    Parameters
    ----------
    point : array_like, shape (2,)
        The query point.
    polygon : array_like, shape (n, 2)
        Polygon vertices (closure point is removed internally).
    include_boundary : bool
        If True, points exactly on the boundary count as inside.
    tol : float
        Numerical tolerance.

    Returns
    -------
    bool
    """
    poly = remove_duplicate_closure(polygon)
    if len(poly) < 3:
        return False

    x, y = point

    # Boundary check first
    for i in range(len(poly)):
        p1 = poly[i]
        p2 = poly[(i + 1) % len(poly)]
        if _point_on_segment_2d(point, p1, p2, tol=tol):
            return include_boundary

    inside = False
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]

        intersects = (y1 > y) != (y2 > y)
        if intersects:
            x_intersect = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x_intersect > x:
                inside = not inside

    return inside


# ---------------------------------------------------------------------------
# Helper: polar projection
# ---------------------------------------------------------------------------
def _project_to_polar_plane(
    coords: np.ndarray,
    pole: PoleType,
) -> np.ndarray:
    """Project lon/lat vertices into a 2-D polar plane.

    The pole maps to the origin; radial distance is ``90 - |lat - pole_lat|``.

    Parameters
    ----------
    coords : numpy.ndarray, shape (n, 2)
        Longitude / latitude array (already validated).
    pole : PoleType
        Which pole to project towards.

    Returns
    -------
    numpy.ndarray, shape (n, 2)
        Projected 2-D coordinates.
    """
    lons = coords[:, 0]
    lats = coords[:, 1]

    lon_rad = np.deg2rad(lons)
    if pole == PoleType.SOUTH:
        radial = np.maximum(0.0, LATITUDE_BOUNDARY + lats)
    else:
        radial = np.maximum(0.0, LATITUDE_BOUNDARY - lats)

    return np.column_stack((radial * np.cos(lon_rad), radial * np.sin(lon_rad)))


# ---------------------------------------------------------------------------
# Module-level convenience function
# ---------------------------------------------------------------------------
def polygon_includes_pole(
    coords: NDArray[np.floating],
    pole: str = "south",
    include_boundary: bool = False,
    tol: float = POLAR_TOLERANCE,
) -> bool:
    """Check whether a polygon includes the requested pole in its interior.

    The test projects lon/lat vertices into a local polar plane where the pole
    maps to the origin, then runs a standard 2-D point-in-polygon query.

    Parameters
    ----------
    coords : array_like, shape (n, 2)
        Polygon vertices as (longitude, latitude) pairs.
    pole : str
        ``"north"`` or ``"south"`` (case-insensitive).
    include_boundary : bool
        If True, a polygon whose boundary passes through the pole counts.
    tol : float
        Numerical tolerance for geometric tests.

    Returns
    -------
    bool
        True if the polygon encloses the specified pole.
    """
    if coords is None:
        return False

    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 3:
        return False

    arr = remove_duplicate_closure(arr)
    if len(arr) < 3:
        return False

    pole_lc = pole.lower()
    try:
        pole_type = PoleType(pole_lc)
    except ValueError:
        raise ValueError("pole must be either 'south' or 'north'")

    projected = _project_to_polar_plane(arr, pole_type)
    origin = np.array([0.0, 0.0])

    return _point_in_polygon_2d(
        origin, projected, include_boundary=include_boundary, tol=tol
    )


# ---------------------------------------------------------------------------
# Main Handler Class
# ---------------------------------------------------------------------------
class PolarHandler:
    """Centralized polar-region detector for polygon geometries.

    Use this class to uniformly test whether polygons enclose a geographic pole.
    Construct with the desired pole and tolerance; then call :meth:`detect` or
    :meth:`includes` on any coordinate array.

    Parameters
    ----------
    pole : str or PoleType
        Which pole to test against.  Accepts ``"north"``, ``"south"``,
        or a :class:`PoleType` enum member.  Default is ``"south"``.
    include_boundary : bool
        If True, a polygon whose boundary passes through the pole counts
        as including it.  Default is False.
    tolerance : float
        Numerical tolerance for geometric tests.  Default is 1e-12.

    Examples
    --------
    >>> import numpy as np
    >>> from pyearth.gis.polar_handler import PolarHandler
    >>> handler = PolarHandler(pole="south")
    >>> coords = np.array([
    ...     [0.0, -80.0], [90.0, -80.0],
    ...     [180.0, -80.0], [-90.0, -80.0],
    ...     [0.0, -80.0],
    ... ])
    >>> result = handler.detect(coords)
    >>> result.includes_pole
    True
    """

    def __init__(
        self,
        pole: str | PoleType = PoleType.SOUTH,
        include_boundary: bool = False,
        tolerance: float = POLAR_TOLERANCE,
    ):
        if isinstance(pole, PoleType):
            self._pole = pole
        else:
            try:
                self._pole = PoleType(pole.lower())
            except (ValueError, AttributeError):
                raise ValueError("pole must be either 'south' or 'north'")

        self.include_boundary = include_boundary
        self.tolerance = tolerance

    # -- public API ---------------------------------------------------------

    def detect(self, coords: NDArray[np.floating]) -> PolarResult:
        """Detect whether a polygon encloses the configured pole.

        Parameters
        ----------
        coords : array_like, shape (n, 2)
            Polygon vertices as (longitude, latitude) pairs.

        Returns
        -------
        PolarResult
            Result object with ``includes_pole``, ``pole``, and
            ``projected_coords`` fields.
        """
        arr = self._validate_coords(coords)
        if arr is None:
            return PolarResult(includes_pole=False, pole=self._pole)

        projected = _project_to_polar_plane(arr, self._pole)
        origin = np.array([0.0, 0.0])

        includes = _point_in_polygon_2d(
            origin,
            projected,
            include_boundary=self.include_boundary,
            tol=self.tolerance,
        )
        return PolarResult(
            includes_pole=includes,
            pole=self._pole,
            projected_coords=projected,
        )

    def includes(self, coords: NDArray[np.floating]) -> bool:
        """Convenience wrapper: return True if the polygon encloses the pole.

        Parameters
        ----------
        coords : array_like, shape (n, 2)
            Polygon vertices as (longitude, latitude) pairs.

        Returns
        -------
        bool
        """
        return self.detect(coords).includes_pole

    # -- private helpers ----------------------------------------------------

    def _validate_coords(self, coords: NDArray[np.floating]) -> Optional[np.ndarray]:
        """Validate and prepare a coordinate array for processing.

        Returns ``None`` when the input is invalid (too few vertices,
        wrong shape, etc.).
        """
        if coords is None:
            return None

        arr = np.asarray(coords, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 3:
            return None

        arr = remove_duplicate_closure(arr)
        if len(arr) < 3:
            return None

        return arr
