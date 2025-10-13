from typing import List, Dict, Any, Tuple
from pyearth.toolbox.mesh.vertex import pyvertex
from pyearth.toolbox.mesh.edge import pyedge
from pyearth.toolbox.mesh.flowline import pyflowline
from pyearth.toolbox.mesh.polygon import pypolygon
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates
from pyearth.gis.geometry.calculate_angle_between_vertex import calculate_angle_between_vertex

def _create_vertex(coords: Dict[str, Any]) -> pyvertex:
    """Helper function to create a pyvertex object."""
    return pyvertex(coords)

def convert_gcs_coordinates_to_cell(iMesh_type_in: int,
                                    dLongitude_center_in: float,
                                    dLatitude_center_in: float,
                                    aCoordinates_gcs_in: List[Tuple[float, float]],
                                    iFlag_simplify_in: bool = False) -> pypolygon:
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
    vertices = [_create_vertex({'dLongitude_degree': lon, 'dLatitude_degree': lat}) for lon, lat in aCoordinates_gcs_in[:-1]]

    if iFlag_simplify_in:
        simplified_vertices = []
        num_vertices = len(vertices)
        for i in range(num_vertices):
            pv_start = vertices[i - 1]
            pv_middle = vertices[i]
            pv_end = vertices[(i + 1) % num_vertices]

            angle3deg = calculate_angle_between_vertex(pv_start.dLongitude_degree, pv_start.dLatitude_degree,
                                                   pv_middle.dLongitude_degree, pv_middle.dLatitude_degree,
                                                   pv_end.dLongitude_degree, pv_end.dLatitude_degree)

            if angle3deg <= 175:
                simplified_vertices.append(pv_middle)
        vertices = simplified_vertices

    edges = [pyedge(vertices[i], vertices[(i + 1) % len(vertices)]) for i in range(len(vertices))]

    return pypolygon(dLongitude_center_in, dLatitude_center_in, edges, vertices)

def convert_gcs_coordinates_to_flowline(aCoordinates_in: List[Tuple[float, float]]) -> pyflowline:
    """
    Convert GCS coordinates to a pyflowline object.

    Args:
        aCoordinates_in (List[Tuple[float, float]]): A list of GCS coordinates.

    Returns:
        pyflowline: A pyflowline object or None if no edge is created.
    """
    vertices = [_create_vertex({'dLongitude_degree': lon, 'dLatitude_degree': lat}) for lon, lat in aCoordinates_in]
    edges = [pyedge(vertices[j], vertices[j + 1]) for j in range(len(vertices) - 1) if vertices[j] != vertices[j + 1]]

    if not edges:
        print('No edge is created')
        return None

    return pyflowline(edges)

def convert_pcs_coordinates_to_flowline(aCoordinates_in: List[Tuple[float, float]], pProjection_in: Any) -> pyflowline:
    """
    Convert PCS coordinates to a pyflowline object.

    Args:
        aCoordinates_in (List[Tuple[float, float]]): A list of PCS coordinates.
        pProjection_in (Any): The projection information.

    Returns:
        pyflowline: A pyflowline object.
    """
    vertices = []
    for x, y in aCoordinates_in:
        lon, lat = reproject_coordinates(x, y, pProjection_in)
        vertices.append(_create_vertex({'x': x, 'y': y, 'lon': lon, 'lat': lat}))

    edges = [pyedge(vertices[j], vertices[j + 1]) for j in range(len(vertices) - 1)]

    return pyflowline(edges)
