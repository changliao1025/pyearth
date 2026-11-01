"""
Clean zero-width spikes, antennae, and retracing boundary lines from OGR geometries.
"""

from typing import List, Optional, Tuple

from osgeo import ogr


def _clean_single_ring(
    ring_geom: Optional[ogr.Geometry], atol: float = 1e-7
) -> Optional[ogr.Geometry]:
    """
    Remove zero-width retracing spikes from a single OGR LinearRing.
    """
    if ring_geom is None or ring_geom.GetPointCount() < 4:
        return ring_geom.Clone() if ring_geom else None

    n_pts = ring_geom.GetPointCount()
    cleaned: List[Tuple[float, float]] = [
        (ring_geom.GetPoint(k)[0], ring_geom.GetPoint(k)[1]) for k in range(n_pts)
    ]
    changed = True

    while changed:
        changed = False
        new_pts: List[Tuple[float, float]] = []
        i = 0
        while i < len(cleaned):
            # Skip consecutive identical points
            if i + 1 < len(cleaned) and (
                abs(cleaned[i][0] - cleaned[i + 1][0]) <= atol
                and abs(cleaned[i][1] - cleaned[i + 1][1]) <= atol
            ):
                changed = True
                i += 1
                continue

            # Check for retracing: new_pts[-1] == cleaned[i+1]
            if (
                len(new_pts) > 0
                and i + 1 < len(cleaned)
                and abs(new_pts[-1][0] - cleaned[i + 1][0]) <= atol
                and abs(new_pts[-1][1] - cleaned[i + 1][1]) <= atol
            ):
                changed = True
                i += 2
                continue

            new_pts.append(cleaned[i])
            i += 1
        cleaned = new_pts

    # Check loop wrap-around at ends
    while (
        len(cleaned) > 3
        and abs(cleaned[1][0] - cleaned[-2][0]) <= atol
        and abs(cleaned[1][1] - cleaned[-2][1]) <= atol
    ):
        cleaned = cleaned[1:-1]
        cleaned.append(cleaned[0])

    if len(cleaned) < 3:
        return ring_geom.Clone()

    if (
        abs(cleaned[0][0] - cleaned[-1][0]) > atol
        or abs(cleaned[0][1] - cleaned[-1][1]) > atol
    ):
        cleaned.append(cleaned[0])

    new_ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in cleaned:
        new_ring.AddPoint(x, y)
    new_ring.CloseRings()
    return new_ring


def clean_ring_spikes(
    geometry: Optional[ogr.Geometry], atol: float = 1e-7
) -> Optional[ogr.Geometry]:
    """
    Remove zero-width spikes, antennae, and retracing boundary segments from OGR geometries.

    Supports LinearRing, Polygon, and MultiPolygon geometry types.

    Parameters
    ----------
    geometry : ogr.Geometry, optional
        Input geometry (LinearRing, Polygon, or MultiPolygon).
    atol : float, default=1e-7
        Absolute tolerance in coordinate units (degrees or meters) for identifying
        identical or retracing vertex positions.

    Returns
    -------
    ogr.Geometry or None
        Cleaned geometry with retracing spikes removed.
    """
    if geometry is None:
        return None

    geom_name = geometry.GetGeometryName().upper()

    if geom_name == "LINEARRING":
        return _clean_single_ring(geometry, atol=atol)

    elif geom_name == "POLYGON":
        poly_geom = ogr.Geometry(ogr.wkbPolygon)
        for i in range(geometry.GetGeometryCount()):
            r = geometry.GetGeometryRef(i)
            clean_r = _clean_single_ring(r, atol=atol)
            if clean_r is not None:
                poly_geom.AddGeometry(clean_r)
        return poly_geom

    elif geom_name == "MULTIPOLYGON":
        multi_geom = ogr.Geometry(ogr.wkbMultiPolygon)
        for i in range(geometry.GetGeometryCount()):
            part = geometry.GetGeometryRef(i)
            clean_part = clean_ring_spikes(part, atol=atol)
            if clean_part is not None:
                multi_geom.AddGeometry(clean_part)
        return multi_geom

    return geometry.Clone()


__all__ = ["clean_ring_spikes"]
