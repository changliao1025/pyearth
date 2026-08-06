"""Unified International Date Line handler for polygon processing.

This module provides a centralized, strategy-based approach to handling polygons
that cross the International Date Line (IDL). It consolidates all IDL-related
operations into a single class with well-defined strategies.

Usage:
    from pyearth.gis.geometry.idl_handler import IdlHandler, IdlStrategy

    handler = IdlHandler(IdlStrategy.SPLIT)

    # Detect crossing
    result = handler.detect(coords_array)

    # Transform OGR geometry
    fixed_geoms = handler.fix(ogr_polygon)

    # Unwrap for area calculations
    unwrapped = handler.unwrap(coords_array)
"""


from __future__ import annotations

import numpy as np
from osgeo import ogr
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple, Union


# ---------------------------------------------------------------------------
# Module-level constants (single source of truth)
# ---------------------------------------------------------------------------
IDL_TOLERANCE = 1e-6   # Tolerance for detecting points on the IDL (+/-180 deg)
IDL_OFFSET = 1e-7      # Offset to move points off IDL (< IDL_TOLERANCE)
DEFAULT_CUTOFF = -150  # Western-hemisphere threshold for longitude shifting


# ---------------------------------------------------------------------------
# Enums & Data classes
# ---------------------------------------------------------------------------
class IdlStrategy(Enum):
    """Strategy for handling polygons that cross the International Date Line.

    Attributes
    ----------
    CHECK_ONLY : str
        Only detect whether the polygon crosses the IDL. No transformation.
    DROP : str
        Discard any polygons that cross the IDL entirely.
    SPLIT : str
        Split crossing polygons into two sub-polygons (eastern + western hemispheres).
    WRAP_TO_360 : str
        Shift western-hemisphere longitudes by +360 so the polygon becomes continuous
        in the [0, 360] range.
    UNWRAP : str
        Offset all longitudes relative to the first vertex so there are no jumps > 180.
        Best used before area/orientation calculations.
    REORDER : str
        Rotate vertex list until a valid (non-self-intersecting) polygon is found.
    """

    CHECK_ONLY = "check_only"
    DROP = "drop"
    SPLIT = "split"
    WRAP_TO_360 = "wrap_360"
    UNWRAP = "unwrap"
    REORDER = "reorder"


@dataclass
class IdlResult:
    """Encapsulates the result of an IDL processing operation.

    Attributes
    ----------
    crosses_idl : bool
        Whether the input geometry crosses the International Date Line.
    coords : numpy.ndarray or None
        Transformed coordinates when strategy involves coordinate manipulation.
    ogr_geometry : ogr.Geometry or None
        Transformed OGR geometry when using SPLIT or WRAP_TO_360.
    sub_polygons : list of numpy.ndarray or None
        For SPLIT strategy: [eastern_polygon, western_polygon].
    adjusted_coords : numpy.ndarray or None
        When IDL-touching vertices were nudged into the dominant hemisphere.
    """

    crosses_idl: bool
    coords: Optional[np.ndarray] = None
    ogr_geometry: Optional[Union[ogr.Geometry, List[ogr.Geometry]]] = None
    sub_polygons: Optional[List[np.ndarray]] = None
    adjusted_coords: Optional[np.ndarray] = None


# ---------------------------------------------------------------------------
# Helper: Shoelace formula (canonical location)
# ---------------------------------------------------------------------------
def calculate_signed_area_shoelace(coords: np.ndarray) -> float:
    """Calculate the signed area of a polygon using the shoelace formula."""
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")
    x, y = coords[:, 0], coords[:, 1]
    x_rolled = np.roll(x, -1)
    y_rolled = np.roll(y, -1)
    return 0.5 * np.sum(x * y_rolled - x_rolled * y)


def _detect_crossing(coords: np.ndarray) -> Tuple[bool, Optional[np.ndarray]]:
    """Internal IDL-crossing detection logic."""
    if not isinstance(coords, np.ndarray) or coords.ndim != 2:
        raise ValueError("coords must be a 2D numpy array")
    if coords.shape[1] != 2:
        raise ValueError("coords must have 2 columns (lon, lat)")
    if len(coords) < 3:
        return False, None
    lons = coords[:, 0]
    idl_vertices = np.abs(np.abs(lons) - 180.0) < IDL_TOLERANCE
    if np.any(idl_vertices):
        lons_next = np.roll(lons, -1)
        eastward_crossings = ((lons > 0) & (lons_next < 0)) & ~idl_vertices & ~np.roll(idl_vertices, -1)
        westward_crossings = ((lons < 0) & (lons_next > 0)) & ~idl_vertices & ~np.roll(idl_vertices, -1)
        touches_both = bool(np.any(lons > 0)) and bool(np.any(lons < 0))
        spans_both = bool(abs(lons.max() - lons.min()) >= 180.0)
        if touches_both or spans_both:
            return True, None
        elif not (np.any(eastward_crossings) or np.any(westward_crossings)):
            coords_updated = coords.copy()
            for idx in np.where(idl_vertices)[0]:
                prev_idx = (idx - 1) % len(coords_updated)
                next_idx = (idx + 1) % len(coords_updated)
                prev_lon = coords_updated[prev_idx, 0]
                next_lon = coords_updated[next_idx, 0]
                neighbor_lons = []
                if abs(abs(prev_lon) - 180.0) > IDL_TOLERANCE:
                    neighbor_lons.append(prev_lon)
                if abs(abs(next_lon) - 180.0) > IDL_TOLERANCE:
                    neighbor_lons.append(next_lon)
                if neighbor_lons:
                    positive_n = sum(1 for lon in neighbor_lons if lon > 0)
                    if positive_n >= len(neighbor_lons) / 2:
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET
                else:
                    if np.sum(lons > 0) >= np.sum(lons < 0):
                        coords_updated[idx, 0] = 180.0 - IDL_OFFSET
                    else:
                        coords_updated[idx, 0] = -180.0 + IDL_OFFSET
            return False, coords_updated
        else:
            return True, None
    else:
        lon_diffs = np.abs(np.diff(lons))
        max_jump = np.max(lon_diffs)
        wrap_jump = abs(lons[-1] - lons[0])
        crossing = bool(max_jump > 180 or wrap_jump > 180)
        return crossing, None


def _unwrap_longitudes(coords: np.ndarray) -> np.ndarray:
    """Unwrap longitudes relative to the first vertex."""
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")
    result = coords.copy()
    lons = result[:, 0]
    ref_lon = lons[0]
    diff = lons - ref_lon
    lons[diff > 180] -= 360
    lons[diff < -180] += 360
    return result


def _shift_western_by_360(coords: np.ndarray, cutoff: float = DEFAULT_CUTOFF) -> np.ndarray:
    """Shift western-hemisphere longitudes by +360 for continuous [0,360] range."""
    result = coords.copy()
    mask = result[:, 0] < cutoff
    result[mask, 0] += 360.0
    return result


def _check_ccw(coords: np.ndarray) -> bool:
    """Check if polygon vertices are counter-clockwise."""
    if not isinstance(coords, np.ndarray) or coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError("coords must be a 2D numpy array with shape (n, 2)")
    if len(coords) < 3:
        return True
    return calculate_signed_area_shoelace(coords) > 0


def _reverse_if_cw(coords: np.ndarray) -> np.ndarray:
    """Reverse vertex order if clockwise, keeping closure point intact."""
    if not _check_ccw(coords):
        unclosed = coords[:-1]
        return np.vstack([unclosed[::-1], coords[-1:]]).copy()
    return coords


def _ensure_valid_range(coords: np.ndarray) -> np.ndarray:
    """Ensure longitudes are in [-180, 180]."""
    result = coords.copy()
    result[result[:, 0] > 180, 0] -= 360
    result[result[:, 0] < -180, 0] += 360
    return result


def _split_international_date_line_polygon(aCoord_gcs: np.ndarray) -> List[np.ndarray]:
    """Split a polygon crossing the IDL into eastern and western sub-polygons."""
    from pyearth.gis.geometry.calculate_intersect_on_great_circle import (
        find_great_circle_intersection_with_meridian,
    )
    coords = np.array(aCoord_gcs, dtype=np.float64)
    if coords.shape[0] < 4:
        raise ValueError("Polygon must have at least 3 unique vertices plus closure")
    if not np.allclose(coords[0], coords[-1]):
        coords = np.vstack([coords, coords[0]])
    if np.any(coords[:, 0] < -180) or np.any(coords[:, 0] > 180):
        raise ValueError("Longitudes must be in [-180, 180]")
    if np.any(coords[:, 1] < -90) or np.any(coords[:, 1] > 90):
        raise ValueError("Latitudes must be in [-90, 90]")
    coords = _reverse_if_cw(coords)
    lons = coords[:, 0]
    lons_next = np.roll(lons, -1)
    idl_vtx = np.abs(np.abs(lons) - 180.0) < IDL_TOLERANCE
    eastward = ((lons > 0) & (lons_next < 0)) & ~idl_vtx & ~np.roll(idl_vtx, -1)
    westward = ((lons < 0) & (lons_next > 0)) & ~idl_vtx & ~np.roll(idl_vtx, -1)
    crossing_edge_indices = (
        np.where(eastward)[0].tolist() + np.where(westward)[0].tolist()
    )
    if len(crossing_edge_indices) != 2:
        if np.any(lons < 0):
            for i in range(len(coords) - 1):
                if abs(abs(lons[i]) - 180.0) < IDL_TOLERANCE:
                    coords[i, 0] = -180.0 + IDL_OFFSET
            western = _ensure_valid_range(coords[:-1].copy())
            eastern = np.empty((0, 2), dtype=np.float64)
        else:
            for i in range(len(coords) - 1):
                if abs(abs(lons[i]) - 180.0) < IDL_TOLERANCE:
                    coords[i, 0] = 180.0 - IDL_OFFSET
            eastern = _ensure_valid_range(coords[:-1].copy())
            western = np.empty((0, 2), dtype=np.float64)
        return [eastern, western]
    lat_map = {}
    for edge_idx in crossing_edge_indices:
        i_cur, i_nxt = edge_idx, (edge_idx + 1) % (len(coords) - 1)
        lon_a, lat_a = coords[i_cur, 0], coords[i_cur, 1]
        lon_b, lat_b = coords[i_nxt, 0], coords[i_nxt, 1]
        intersection_result = find_great_circle_intersection_with_meridian(
            lon_a, lat_a, lon_b, lat_b
        )
        if intersection_result is not None:
            lat_map[edge_idx] = intersection_result
    EASTERN_BOUNDARY = 180.0 - IDL_OFFSET
    WESTERN_BOUNDARY = -180.0 + IDL_OFFSET
    def _boundary_for(sub_list):
        if not sub_list:
            return WESTERN_BOUNDARY
        if any(pt[0] > 0 for pt in sub_list):
            return EASTERN_BOUNDARY
        return WESTERN_BOUNDARY
    def _snap_append(lst, lon, lat):
        if abs(abs(lon) - 180.0) < IDL_TOLERANCE:
            lon = _boundary_for(lst)
        lst.append([lon, lat])
    eastern_pts: list = []
    western_pts: list = []
    active = eastern_pts
    inactive = western_pts
    for i in range(len(coords) - 1):
        _snap_append(active, coords[i, 0], coords[i, 1])
        if i in lat_map:
            lat_cross = lat_map[i]
            active.append([_boundary_for(active), lat_cross])
            active, inactive = inactive, active
            active.append([_boundary_for(active), lat_cross])
    for pts in (eastern_pts, western_pts):
        if len(pts) >= 2 and not np.allclose(pts[0], pts[-1]):
            pts.append(list(pts[0]))
    eastern_arr = _ensure_valid_range(np.array(eastern_pts)) if eastern_pts else np.empty((0, 2), dtype=np.float64)
    western_arr = _ensure_valid_range(np.array(western_pts)) if western_pts else np.empty((0, 2), dtype=np.float64)
    if eastern_arr.size:
        eastern_arr = _reverse_if_cw(eastern_arr)
    if western_arr.size:
        western_arr = _reverse_if_cw(western_arr)
    return [eastern_arr, western_arr]


def _reorder_vertices_until_valid(vertices: list) -> list:
    """Reorder polygon vertices by rotation until OGR reports validity."""
    if not vertices:
        return []
    for i, v in enumerate(vertices):
        lon, lat = float(v[0]), float(v[1])
        if not (-180 <= lon <= 180):
            raise ValueError(f"Longitude {lon} at vertex {i} out of range")
        if not (-90 <= lat <= 90):
            raise ValueError(f"Latitude {lat} at vertex {i} out of range")
    if len(vertices) < 3:
        raise ValueError("Polygon must have at least 3 vertices")
    if vertices[0] != vertices[-1]:
        vertices = vertices + [vertices[0]]
    n_points = len(vertices) - 1
    current = vertices.copy()
    for _ in range(n_points + 1):
        try:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for lon, lat in current:
                ring.AddPoint(float(lon), float(lat))
            ring.CloseRings()
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            if poly.IsValid():
                return current + [current[0]]
        except Exception:
            pass
        current = current[1:] + [current[0]]
    raise ValueError(
        "Cannot form a valid polygon by rotating vertices. "
        "The input may be self-intersecting or invalid."
    )


def _convert_to_unwrapped_polygon(geometry_in: ogr.Geometry) -> Optional[ogr.Geometry]:
    """Convert an IDL-crossing OGR polygon by shifting western longitudes +360."""
    if geometry_in is None:
        raise ValueError("Input geometry cannot be None")
    srs = None
    try:
        srs = geometry_in.CloneSRG()
    except Exception:
        srs = None
    n_geom = geometry_in.GetGeometryCount()
    if n_geom > 0 and geometry_in.GetGeometryType() == ogr.wkbMultiPolygon:
        parts = []
        for i in range(n_geom):
            part = geometry_in.GetGeometryRef(i)
            if part is not None:
                converted = _convert_single_idl_polygon(part)
                if converted is not None:
                    parts.append(converted)
        if parts:
            multi = ogr.Geometry(ogr.wkbMultiPolygon)
            for p in parts:
                multi.AddGeometry(p)
            return multi if len(parts) > 1 else parts[0]
        return None
    return _convert_single_idl_polygon(geometry_in)


def _convert_single_idl_polygon(geometry_in: ogr.Geometry) -> Optional[ogr.Geometry]:
    """Convert a single IDL-crossing polygon by shifting western longitudes +360."""
    if geometry_in.GetGeometryType() != ogr.wkbPolygon:
        return None
    outer = geometry_in.GetGeometryRef(0)
    if outer is None:
        return None
    n_pts = outer.GetPointCount()
    coords = [(outer.GetX(i), outer.GetY(i)) for i in range(n_pts)]
    lons = [c[0] for c in coords]
    span = max(lons) - min(lons)
    if span < 180:
        if geometry_in.IsValid():
            return geometry_in
        return None
    new_coords = []
    for lon, lat in coords:
        if lon < DEFAULT_CUTOFF:
            lon += 360.0
        new_coords.append((lon, lat))
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for lon, lat in new_coords:
        ring.AddPoint(lon, lat)
    ring.CloseRings()
    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)
    if poly.IsValid():
        if srs is not None:
            poly.AssignSpatialReference(srs)
        return poly
    return None


# ---------------------------------------------------------------------------
# Main Handler Class
# ---------------------------------------------------------------------------
class IdlHandler:
    """Centralized International Date Line processor for polygon geometries.

    Use this class to uniformly handle IDL-crossing polygons across your codebase.
    Choose a strategy at construction time; all subsequent operations follow that
    strategy unless overridden.

    Parameters
    ----------
    strategy : IdlStrategy
        How to handle IDL-crossing polygons. Default is ``CHECK_ONLY``.
    longitude_cutoff : float
        The threshold below which longitudes are considered "western".
        Default is -150.
    idl_tolerance : float
        Numerical tolerance for +/-180 detection. Default is 1e-6.
    idl_offset : float
        Offset applied when nudging vertices off +/-180. Default is 1e-7.

    Examples
    --------
    >>> from pyearth.gis.geometry.idl_handler import IdlHandler, IdlStrategy
    >>> handler = IdlHandler(IdlStrategy.SPLIT)
    >>> result = handler.detect(np.array([[170, -10], [170, 10], [-170, 10], [-170, -10], [170, -10]]))
    >>> result.crosses_idl
    True
    """

    def __init__(
        self,
        strategy: IdlStrategy = IdlStrategy.CHECK_ONLY,
        longitude_cutoff: float = DEFAULT_CUTOFF,
        idl_tolerance: float = IDL_TOLERANCE,
        idl_offset: float = IDL_OFFSET,
    ):
        self.strategy = strategy
        self.longitude_cutoff = longitude_cutoff
        self.idl_tolerance = idl_tolerance
        self.idl_offset = idl_offset

    def detect(self, coords: np.ndarray) -> IdlResult:
        """Detect whether a polygon crosses the IDL."""
        crosses, adjusted = _detect_crossing(coords)
        result = IdlResult(crosses_idl=crosses, adjusted_coords=adjusted)
        if crosses and self.strategy == IdlStrategy.SPLIT:
            result.sub_polygons = _split_international_date_line_polygon(coords)
        return result

    def fix(self, ogr_geom: ogr.Geometry) -> List[ogr.Geometry]:
        """Fix an OGR geometry that may cross the IDL based on configured strategy."""
        if ogr_geom is None:
            return []
        try:
            coords = self._extract_coords(ogr_geom)
        except (ValueError, AttributeError):
            return [ogr_geom]
        crosses, _ = _detect_crossing(coords)
        if not crosses:
            return [ogr_geom]
        if self.strategy == IdlStrategy.DROP:
            return []
        if self.strategy == IdlStrategy.SPLIT:
            sub_polys = _split_international_date_line_polygon(coords)
            results = []
            for sp in sub_polys:
                if sp.size > 0 and len(sp) >= 3:
                    results.append(self._coords_to_ogr(sp, ogr_geom.CloneSRG()))
            return results
        if self.strategy == IdlStrategy.WRAP_TO_360:
            unwrapped = _shift_western_by_360(coords, self.longitude_cutoff)
            return [self._coords_to_ogr(unwrapped, ogr_geom.CloneSRG())]
        return [ogr_geom]

    def unwrap(self, coords: np.ndarray) -> np.ndarray:
        """Unwrap longitudes so they are continuous (for area calculations)."""
        return _unwrap_longitudes(coords)

    @staticmethod
    def check_counter_clockwise(coords: np.ndarray) -> bool:
        """Check if polygon coordinates are in counter-clockwise order.

        Handles IDL crossings by unwrapping longitudes first if needed.
        """
        crosses, _ = _detect_crossing(coords)
        if crosses:
            unwrapped = _unwrap_longitudes(coords)
            return _check_ccw(unwrapped)
        return _check_ccw(coords)

    @staticmethod
    def calculate_signed_area(coords: np.ndarray) -> float:
        """Calculate the signed area using the shoelace formula."""
        return calculate_signed_area_shoelace(coords)

    @staticmethod
    def _extract_coords(geom: ogr.Geometry) -> np.ndarray:
        """Extract (lon, lat) numpy array from an OGR geometry."""
        from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
        coords = get_geometry_coordinates(geom)
        if coords is None or len(coords) < 3:
            raise ValueError("Insufficient coordinates extracted from geometry")
        return coords

    @staticmethod
    def _coords_to_ogr(coords: np.ndarray, srs=None) -> ogr.Geometry:
        """Convert a numpy (lon, lat) array to an OGR Polygon."""
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for lon, lat in coords[:-1]:
            ring.AddPoint(lon, lat)
        ring.CloseRings()
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        if srs is not None:
            poly.AssignSpatialReference(srs)
        return poly

