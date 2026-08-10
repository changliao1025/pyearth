import numpy as np
from numpy.typing import NDArray
from pyearth.toolbox.data.remove_duplicate_closure import remove_duplicate_closure as remove_duplicate_closure

def _point_on_segment_2d(
    point: NDArray[np.floating],
    seg_start: NDArray[np.floating],
    seg_end: NDArray[np.floating],
    tol: float = 1.0e-12,
) -> bool:
    """Return True if a 2D point lies on a line segment within tolerance."""
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
    tol: float = 1.0e-12,
) -> bool:
    """2D point-in-polygon using ray casting."""
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

def polygon_includes_pole(
    coords: NDArray[np.floating],
    pole: str = "south",
    include_boundary: bool = False,
    tol: float = 1.0e-12,
) -> bool:
    """Check whether a polygon includes the requested pole in its interior.

    The test projects lon/lat vertices into a local polar plane where the pole
    maps to the origin, then runs a standard 2D point-in-polygon query.
    """
    if coords is None:
        return False

    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2 or len(arr) < 3:
        return False

    arr = remove_duplicate_closure(arr)
    if len(arr) < 3:
        return False

    lons = arr[:, 0]
    lats = arr[:, 1]

    lon_rad = np.deg2rad(lons)
    pole_lc = pole.lower()
    if pole_lc == "south":
        radial = np.maximum(0.0, 90.0 + lats)
    elif pole_lc == "north":
        radial = np.maximum(0.0, 90.0 - lats)
    else:
        raise ValueError("pole must be either 'south' or 'north'")

    projected = np.column_stack((radial * np.cos(lon_rad), radial * np.sin(lon_rad)))
    origin = np.array([0.0, 0.0])

    return _point_in_polygon_2d(
        origin, projected, include_boundary=include_boundary, tol=tol
    )