"""
Extract coordinate arrays from OGR geometry objects.

This module provides utilities for extracting coordinate arrays from various
OGR geometry types. For polygon geometries, coordinates are automatically
enforced to be in counter-clockwise (CCW) order following the right-hand rule.
"""

from typing import List, Union
import numpy as np
from osgeo import ogr
from pyearth.gis.geometry.check_counter_clockwise import check_counter_clockwise


def get_geometry_coordinates(
    geometry: ogr.Geometry,
) -> Union[np.ndarray, List[np.ndarray]]:
    """Extract coordinates from an OGR geometry object.

    Dispatches to appropriate handler based on geometry type. For polygon
    geometries, ensures coordinates are in counter-clockwise (CCW) order.

    Parameters
    ----------
    geometry : ogr.Geometry
        OGR geometry object from which to extract coordinates.

    Returns
    -------
    np.ndarray or List[np.ndarray]
        For POINT, LINESTRING, POLYGON, LINEARRING: numpy array of shape (n, 2)
        For MULTIPOLYGON: list of numpy arrays, one per polygon part

    Raises
    ------
    ValueError
        If geometry is None, invalid, or of unsupported type.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    Supported geometry types:
    - POINT: Returns single point as array of shape (1, 2)
    - LINESTRING: Returns array of line vertices
    - LINEARRING: Returns array of ring vertices
    - POLYGON: Returns exterior ring coordinates in CCW order
    - MULTIPOLYGON: Returns list of exterior rings, each in CCW order

    Currently unsupported (will raise ValueError):
    - MULTIPOINT, MULTILINESTRING, GEOMETRYCOLLECTION

    Z-coordinates (3D) are ignored; only X and Y are extracted.

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a point geometry
    >>> point = ogr.Geometry(ogr.wkbPoint)
    >>> point.AddPoint(10.0, 20.0)
    >>> coords = get_geometry_coordinates(point)
    >>> coords.shape
    (1, 2)

    See Also
    --------
    get_polygon_exterior_coords : Extract polygon exterior ring
    get_multipolygon_exterior_coords : Extract multipolygon exterior rings
    """
    # Validate geometry
    if geometry is None:
        raise ValueError("Geometry object cannot be None")

    try:
        geometry_type = geometry.GetGeometryName()
    except Exception as e:
        raise ValueError(f"Invalid geometry object: {e}")

    # Dispatch to appropriate handler
    if geometry_type == "POINT":
        return get_point_coords(geometry)
    elif geometry_type == "LINESTRING":
        return get_linestring_coords(geometry)
    elif geometry_type == "POLYGON":
        return get_polygon_exterior_coords(geometry)
    elif geometry_type == "LINEARRING":
        return get_linearring_coords(geometry)
    elif geometry_type == "MULTIPOLYGON":
        return get_multipolygon_exterior_coords(geometry)
    else:
        raise ValueError(
            f"Unsupported geometry type: '{geometry_type}'. "
            "Supported types: POINT, LINESTRING, POLYGON, LINEARRING, MULTIPOLYGON"
        )


def get_polygon_exterior_coords(polygon_geometry: ogr.Geometry) -> np.ndarray:
    """Extract exterior ring coordinates from a polygon, enforcing CCW order.

    Parameters
    ----------
    polygon_geometry : ogr.Geometry
        OGR Polygon geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing exterior ring coordinates in
        counter-clockwise order.

    Raises
    ------
    ValueError
        If geometry is not a valid polygon or has no exterior ring.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    - Only extracts the exterior ring; holes/interior rings are ignored
    - Automatically reverses coordinate order if clockwise
    - Returns (longitude, latitude) or (x, y) coordinate pairs

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a square polygon
    >>> ring = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring.AddPoint(0, 0)
    >>> ring.AddPoint(1, 0)
    >>> ring.AddPoint(1, 1)
    >>> ring.AddPoint(0, 1)
    >>> ring.AddPoint(0, 0)
    >>> polygon = ogr.Geometry(ogr.wkbPolygon)
    >>> polygon.AddGeometry(ring)
    >>> coords = get_polygon_exterior_coords(polygon)
    >>> coords.shape[1]
    2
    """
    if polygon_geometry is None:
        raise ValueError("Polygon geometry cannot be None")

    try:
        # Get the exterior ring (index 0)
        ring = polygon_geometry.GetGeometryRef(0)
        if ring is None:
            raise ValueError("Polygon has no exterior ring")

        n_points = ring.GetPointCount()
        if n_points == 0:
            raise ValueError("Polygon exterior ring has no points")

        # Extract coordinates
        exterior_coords = []
        for i in range(n_points):
            point = ring.GetPoint(i)
            exterior_coords.append((point[0], point[1]))

    except Exception as e:
        raise RuntimeError(f"Failed to extract polygon coordinates: {e}")

    # Convert to numpy array
    coords_array = np.array(exterior_coords)

    # Ensure counter-clockwise order
    if not check_counter_clockwise(coords_array):
        coords_array = coords_array[::-1]

    return coords_array


def get_multipolygon_exterior_coords(
    multipolygon_geometry: ogr.Geometry,
) -> List[np.ndarray]:
    """Extract exterior ring coordinates from all parts of a multipolygon.

    Each polygon part's coordinates are enforced to be in counter-clockwise order.

    Parameters
    ----------
    multipolygon_geometry : ogr.Geometry
        OGR MultiPolygon geometry object.

    Returns
    -------
    List[np.ndarray]
        List of numpy arrays, one per polygon part. Each array has shape (n, 2)
        and contains coordinates in counter-clockwise order.

    Raises
    ------
    ValueError
        If geometry is not a valid multipolygon or has no parts.
    RuntimeError
        If coordinate extraction fails for any part.

    Notes
    -----
    - Only extracts exterior rings; holes/interior rings are ignored
    - Each polygon part is independently enforced to CCW order
    - Empty parts are skipped with a warning

    Examples
    --------
    >>> from osgeo import ogr
    >>> # Create a multipolygon with 2 parts
    >>> multi = ogr.Geometry(ogr.wkbMultiPolygon)
    >>> # Add first polygon
    >>> ring1 = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring1.AddPoint(0, 0)
    >>> ring1.AddPoint(1, 0)
    >>> ring1.AddPoint(1, 1)
    >>> ring1.AddPoint(0, 0)
    >>> poly1 = ogr.Geometry(ogr.wkbPolygon)
    >>> poly1.AddGeometry(ring1)
    >>> multi.AddGeometry(poly1)
    >>> coords_list = get_multipolygon_exterior_coords(multi)
    >>> len(coords_list)
    1
    """
    if multipolygon_geometry is None:
        raise ValueError("MultiPolygon geometry cannot be None")

    # Validate geometry type
    try:
        geometry_type = multipolygon_geometry.GetGeometryName()
        if geometry_type != "MULTIPOLYGON":
            raise ValueError(f"Expected MULTIPOLYGON geometry, got '{geometry_type}'")
    except Exception as e:
        raise ValueError(f"Invalid geometry object: {e}")

    try:
        n_parts = multipolygon_geometry.GetGeometryCount()
        if n_parts == 0:
            raise ValueError("MultiPolygon has no parts")

        exterior_coords_list = []

        for i in range(n_parts):
            # Get polygon part
            polygon = multipolygon_geometry.GetGeometryRef(i)
            if polygon is None:
                continue  # Skip invalid parts

            # Get exterior ring
            ring = polygon.GetGeometryRef(0)
            if ring is None:
                continue  # Skip polygons without rings

            n_points = ring.GetPointCount()
            if n_points == 0:
                continue  # Skip empty rings

            # Extract coordinates for this polygon part
            part_coords = []
            for j in range(n_points):
                point = ring.GetPoint(j)
                part_coords.append((point[0], point[1]))

            # Convert to numpy array
            coords_array = np.array(part_coords)

            # Ensure counter-clockwise order
            if not check_counter_clockwise(coords_array):
                coords_array = coords_array[::-1]

            exterior_coords_list.append(coords_array)

    except Exception as e:
        raise RuntimeError(f"Failed to extract multipolygon coordinates: {e}")

    if len(exterior_coords_list) == 0:
        raise ValueError("No valid polygon parts found in multipolygon")

    return exterior_coords_list


def get_linestring_coords(linestring_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a linestring geometry.

    Parameters
    ----------
    linestring_geometry : ogr.Geometry
        OGR LineString geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing linestring vertex coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid linestring or has no points.
    RuntimeError
        If coordinate extraction fails.

    Examples
    --------
    >>> from osgeo import ogr
    >>> line = ogr.Geometry(ogr.wkbLineString)
    >>> line.AddPoint(0, 0)
    >>> line.AddPoint(1, 1)
    >>> coords = get_linestring_coords(line)
    >>> coords.shape
    (2, 2)
    """
    if linestring_geometry is None:
        raise ValueError("LineString geometry cannot be None")

    try:
        n_points = linestring_geometry.GetPointCount()
        if n_points == 0:
            raise ValueError("LineString has no points")

        # Extract coordinates using list comprehension
        coords = [
            (linestring_geometry.GetPoint(i)[0], linestring_geometry.GetPoint(i)[1])
            for i in range(n_points)
        ]

    except Exception as e:
        raise RuntimeError(f"Failed to extract linestring coordinates: {e}")

    return np.array(coords)


def get_point_coords(point_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a point geometry.

    Parameters
    ----------
    point_geometry : ogr.Geometry
        OGR Point geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (1, 2) containing the point coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid point.
    RuntimeError
        If coordinate extraction fails.

    Examples
    --------
    >>> from osgeo import ogr
    >>> point = ogr.Geometry(ogr.wkbPoint)
    >>> point.AddPoint(10.5, 20.3)
    >>> coords = get_point_coords(point)
    >>> coords.shape
    (1, 2)
    >>> coords[0]
    array([10.5, 20.3])
    """
    if point_geometry is None:
        raise ValueError("Point geometry cannot be None")

    try:
        point = point_geometry.GetPoint()
        if point is None or len(point) < 2:
            raise ValueError("Invalid point geometry")

    except Exception as e:
        raise RuntimeError(f"Failed to extract point coordinates: {e}")

    return np.array([(point[0], point[1])])


def get_linearring_coords(linearring_geometry: ogr.Geometry) -> np.ndarray:
    """Extract coordinates from a linear ring geometry.

    Parameters
    ----------
    linearring_geometry : ogr.Geometry
        OGR LinearRing geometry object.

    Returns
    -------
    np.ndarray
        Array of shape (n, 2) containing ring vertex coordinates.

    Raises
    ------
    ValueError
        If geometry is not a valid linear ring or has no points.
    RuntimeError
        If coordinate extraction fails.

    Notes
    -----
    - A linear ring is a closed linestring (first point = last point)
    - This function does not enforce CCW order (use for non-polygon contexts)

    Examples
    --------
    >>> from osgeo import ogr
    >>> ring = ogr.Geometry(ogr.wkbLinearRing)
    >>> ring.AddPoint(0, 0)
    >>> ring.AddPoint(1, 0)
    >>> ring.AddPoint(1, 1)
    >>> ring.AddPoint(0, 0)
    >>> coords = get_linearring_coords(ring)
    >>> coords.shape
    (4, 2)
    """
    if linearring_geometry is None:
        raise ValueError("LinearRing geometry cannot be None")

    try:
        n_points = linearring_geometry.GetPointCount()
        if n_points == 0:
            raise ValueError("LinearRing has no points")

        # Extract coordinates using list comprehension
        coords = [
            (linearring_geometry.GetPoint(i)[0], linearring_geometry.GetPoint(i)[1])
            for i in range(n_points)
        ]

    except Exception as e:
        raise RuntimeError(f"Failed to extract linear ring coordinates: {e}")

    return np.array(coords)
