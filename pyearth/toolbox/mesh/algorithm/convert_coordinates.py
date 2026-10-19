from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from pyearth.gis.geometry.calculate_angle_between_point import (
    calculate_angle_between_point,
)
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
from pyearth.toolbox.mesh.point import pypoint

if TYPE_CHECKING:
    from pyearth.toolbox.mesh.line import pyline
    from pyearth.toolbox.mesh.polygon import pypolygon
    from pyearth.toolbox.mesh.polyline import pypolyline


def _create_point(coords: dict[str, Any]) -> pypoint:
    """Helper function to create a pypoint object."""
    return pypoint(coords)


def convert_gcs_coordinates_to_meshcell(
    iMesh_type_in: int,
    dLongitude_center_in: float,
    dLatitude_center_in: float,
    aCoordinates_gcs_in: list[tuple[float, float]],
    iFlag_simplify_in: bool = False,
) -> pypolygon:
    """
    Convert GCS coordinates to a pypolygon object.

    Args:
        iMesh_type_in (int): The mesh type.
        dLongitude_center_in (float): The longitude of the center.
        dLatitude_center_in (float): The latitude of the center.
        aCoordinates_gcs_in (List[Tuple[float, float]]): A list of GCS coordinates.
        iFlag_simplify_in (bool, optional): A flag to simplify the polygon. Defaults to False.

    Returns:
        pypolygon: A pypolygon object.
    """
    from pyearth.toolbox.mesh.line import pyline
    from pyearth.toolbox.mesh.polygon import pypolygon

    points = [
        _create_point({"dLongitude_degree": lon, "dLatitude_degree": lat})
        for lon, lat in aCoordinates_gcs_in[:-1]
    ]

    if iFlag_simplify_in:
        simplified_points = []
        num_points = len(points)
        for i in range(num_points):
            pt_start = points[i - 1]
            pt_middle = points[i]
            pt_end = points[(i + 1) % num_points]

            angle3deg = calculate_angle_between_point(
                pt_start.dLongitude_degree,
                pt_start.dLatitude_degree,
                pt_middle.dLongitude_degree,
                pt_middle.dLatitude_degree,
                pt_end.dLongitude_degree,
                pt_end.dLatitude_degree,
            )

            if angle3deg <= 175:
                simplified_points.append(pt_middle)
        points = simplified_points

    edges = [
        pyline(points[i], points[(i + 1) % len(points)]) for i in range(len(points))
    ]

    # add the closing point back to the points list
    points.append(points[0])

    return pypolygon(dLongitude_center_in, dLatitude_center_in, edges, points)


def convert_gcs_coordinates_to_polyline(
    aCoordinates_in: list[tuple[float, float]],
) -> pypolyline | None:
    """
    Convert GCS coordinates to a pypolyline object.

    Args:
        aCoordinates_in (List[Tuple[float, float]]): A list of GCS coordinates.

    Returns:
        pypolyline: A pypolyline object or None if no edge is created.
    """
    from pyearth.toolbox.mesh.line import pyline
    from pyearth.toolbox.mesh.polyline import pypolyline

    points = [
        _create_point({"dLongitude_degree": lon, "dLatitude_degree": lat})
        for lon, lat in aCoordinates_in
    ]
    edges = [
        pyline(points[j], points[j + 1])
        for j in range(len(points) - 1)
        if points[j] != points[j + 1]
    ]

    if not edges:
        return None

    return pypolyline(edges)


def convert_pcs_coordinates_to_polyline(
    aCoordinates_in: list[tuple[float, float]], pProjection_in: Any
) -> pypolyline:
    """
    Convert PCS coordinates to a pypolyline object.

    Args:
        aCoordinates_in (List[Tuple[float, float]]): A list of PCS coordinates.
        pProjection_in (Any): The projection information.

    Returns:
        pypolyline: A pypolyline object.
    """
    from pyearth.toolbox.mesh.line import pyline
    from pyearth.toolbox.mesh.polyline import pypolyline

    points = []
    for x, y in aCoordinates_in:
        lon, lat = reproject_coordinates(x, y, pProjection_in)
        points.append(_create_point({"x": x, "y": y, "lon": lon, "lat": lat}))

    edges = [pyline(points[j], points[j + 1]) for j in range(len(points) - 1)]

    return pypolyline(edges)
