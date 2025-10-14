import os
import json
from json import JSONEncoder
import numpy as np
from typing import List, Tuple

from pyearth.toolbox.mesh.vertex import pyvertex
from pyearth.toolbox.mesh.circle import pycircle

from pyearth.gis.geometry.calculate_intersect_on_great_circle import calculate_intersect_on_great_circle
from pyearth.gis.gdal.write.vector.gdal_export_vertex_to_vector_file import export_vertex_as_polygon
from pyearth.toolbox.mesh.algorithm.split_by_length import split_edge_by_length
from pyearth.toolbox.mesh.algorithm.find_minimal_enclosing_polygon import find_minimal_enclosing_polygon
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude
from pyearth.gis.geometry.calculate_angle_between_vertex import calculate_angle_between_vertex
from pyearth.gis.geometry.calculate_distance_to_plane import calculate_distance_to_plane

class EdgeClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return super().default(obj)

class pyedge(object):
    """The pyedge class represents an edge between two vertices."""

    def __init__(self, pVertex_start_in: pyvertex, pVertex_end_in: pyvertex):
        """
        Initialize a pyedge object.

        Args:
            pVertex_start_in (pyvertex): The starting vertex.
            pVertex_end_in (pyvertex): The ending vertex.

        Raises:
            ValueError: If the two vertices are the same.
        """
        if pVertex_start_in == pVertex_end_in:
            raise ValueError("The two vertices are the same.")

        self.pVertex_start: pyvertex = pVertex_start_in
        self.pVertex_end: pyvertex = pVertex_end_in
        self.dLength: float = self.calculate_length()
        self.pBound: Tuple[float, float, float, float] = self.calculate_edge_bound()
        self.lEdgeID: int = -1
        self.lEdgeIndex: int = -1
        self.lIndex_upstream: int = -1
        self.lIndex_downstream: int = -1

    @classmethod
    def create(cls, pVertex_start_in: pyvertex, pVertex_end_in: pyvertex) -> "pyedge":
        """
        Factory method to create a pyedge object.

        Args:
            pVertex_start_in (pyvertex): The starting vertex.
            pVertex_end_in (pyvertex): The ending vertex.

        Returns:
            pyedge: A pyedge object.
        """
        return cls(pVertex_start_in, pVertex_end_in)

    def calculate_edge_bound(self) -> Tuple[float, float, float, float]:
        """Calculate the bounding box of the edge."""
        dLon_max = max(self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree)
        dLon_min = min(self.pVertex_start.dLongitude_degree, self.pVertex_end.dLongitude_degree)
        dLat_max = max(self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree)
        dLat_min = min(self.pVertex_start.dLatitude_degree, self.pVertex_end.dLatitude_degree)
        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def calculate_length(self) -> float:
        """Calculate the length of the edge."""
        self.dLength = self.pVertex_start.calculate_distance(self.pVertex_end)
        return self.dLength

    def check_shared_vertex(self, other: "pyedge") -> bool:
        """Check whether two edges are sharing the same vertex."""
        return self.pVertex_start in (other.pVertex_start, other.pVertex_end) or \
               self.pVertex_end in (other.pVertex_start, other.pVertex_end)

    def check_upstream(self, other: "pyedge") -> bool:
        """Check whether another edge is the upstream of current edge."""
        return self.pVertex_start == other.pVertex_end

    def check_downstream(self, other: "pyedge") -> bool:
        """Check whether another edge is the downstream of current edge."""
        return self.pVertex_end == other.pVertex_start

    def split_by_length(self, dLength_in: float) -> List["pyedge"]:
        """Split an edge using the threshold."""
        if self.dLength <= dLength_in:
            return [self]
        else:
            return split_edge_by_length(self, dLength_in)

    def reverse(self):
        """Reverse an edge."""
        self.pVertex_start, self.pVertex_end = self.pVertex_end, self.pVertex_start

    def is_overlap(self, pEdge_in: "pyedge") -> bool:
        """Check if two edges overlap each other."""
        return (self.pVertex_start == pEdge_in.pVertex_start and self.pVertex_end == pEdge_in.pVertex_end) or \
               (self.pVertex_start == pEdge_in.pVertex_end and self.pVertex_end == pEdge_in.pVertex_start)

    def check_vertex_on_edge(self, pVertex_in: pyvertex) -> Tuple[bool, float, float]:
        """Check if a vertex on an edge."""
        dDistance = -1.0
        dDistance_plane = 9999.0
        on_edge = False

        if pVertex_in != self.pVertex_start and pVertex_in != self.pVertex_end:
            d1 = self.pVertex_start.calculate_distance(pVertex_in)
            d2 = self.pVertex_end.calculate_distance(pVertex_in)
            d3 = d1 + d2 - self.dLength
            angle3deg = calculate_angle_between_vertex(
                self.pVertex_start.dLongitude_degree, self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree, self.pVertex_end.dLatitude_degree)

            dDistance_plane = calculate_distance_to_plane(
                self.pVertex_start.dLongitude_degree, self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree, self.pVertex_end.dLatitude_degree)

            if angle3deg > 178 and d3 < 1.0:
                on_edge = True
                dDistance = d1

        return on_edge, dDistance, dDistance_plane

    def calculate_distance_to_vertex(self, pVertex_in: pyvertex) -> Tuple[float, pyvertex]:
        """Calculate the minimum distance from a vertex to the edge."""
        d1 = self.pVertex_start.calculate_distance(pVertex_in)
        d2 = self.pVertex_end.calculate_distance(pVertex_in)

        angle3deg = calculate_angle_between_vertex(
            self.pVertex_start.dLongitude_degree, self.pVertex_start.dLatitude_degree,
            pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,
            self.pVertex_end.dLongitude_degree, self.pVertex_end.dLatitude_degree)

        if angle3deg > 90:
            dDistance_plane = calculate_distance_to_plane(
                self.pVertex_start.dLongitude_degree, self.pVertex_start.dLatitude_degree,
                pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,
                self.pVertex_end.dLongitude_degree, self.pVertex_end.dLatitude_degree)

            if dDistance_plane < d1 and dDistance_plane < d2:
                dLongitude_intersect, dLatitude_intersect = calculate_intersect_on_great_circle(
                    self.pVertex_start.dLongitude_degree, self.pVertex_start.dLatitude_degree,
                    pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,
                    self.pVertex_end.dLongitude_degree, self.pVertex_end.dLatitude_degree)

                pVertex_out = pyvertex({'dLongitude_degree': dLongitude_intersect, 'dLatitude_degree': dLatitude_intersect})
                dDistance_min = pVertex_out.calculate_distance(pVertex_in)

                if dDistance_min < d1 and dDistance_min < d2:
                    return dDistance_min, pVertex_out

        if d1 < d2:
            return d1, self.pVertex_start
        else:
            return d2, self.pVertex_end

    def calculate_buffer_zone_polygon(self, dRadius: float, nPoint: int = 36, sFilename_out: str = None, sFolder_out: str = None) -> Tuple[List[pyvertex], List[pyvertex], List[pyvertex], List[pycircle]]:
        """Calculate the buffer zone polygon of the edge."""
        if self.dLength < dRadius * 2.0:
            aEdge = [self]
        else:
            aEdge = self.split_by_length(dRadius)

        aVertex_out = []
        aVertex_center = []
        aVertex_circle = []
        aCircle = []

        for i, pEdge in enumerate(aEdge):
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end

            aVertex_center.extend([pVertex_start, pVertex_end])

            aVertex_start_buffer = pVertex_start.calculate_buffer_zone_circle(dRadius, nPoint)
            pEdge.pCircle_start = pycircle(pVertex_start, aVertex_start_buffer)
            aVertex_circle.extend(aVertex_start_buffer)

            aVertex_end_buffer = pVertex_end.calculate_buffer_zone_circle(dRadius, nPoint)
            pEdge.pCircle_end = pycircle(pVertex_end, aVertex_end_buffer)
            aVertex_circle.extend(aVertex_end_buffer)

            if sFolder_out:
                export_vertex_as_polygon(aVertex_start_buffer, os.path.join(sFolder_out, f'buffer_zone_start_{i}.geojson'))
                export_vertex_as_polygon(aVertex_end_buffer, os.path.join(sFolder_out, f'buffer_zone_end_{i}.geojson'))

            aCircle.extend([pEdge.pCircle_start, pEdge.pCircle_end])

        aLongitude_degree = [v.dLongitude_degree for v in aVertex_circle]
        aLatitude_degree = [v.dLatitude_degree for v in aVertex_circle]

        pPolygon_out = find_minimal_enclosing_polygon(aLongitude_degree, aLatitude_degree)

        aVertex_out = [pyvertex({'dLongitude_degree': p[0], 'dLatitude_degree': p[1]}) for p in pPolygon_out]

        if sFilename_out:
            export_vertex_as_polygon(aVertex_out, sFilename_out)

        return aVertex_out, aVertex_center, aVertex_circle, aCircle

    def __eq__(self, other: object) -> bool:
        """Check if two edges are equivalent."""
        if not isinstance(other, pyedge):
            return NotImplemented
        return self.pVertex_start == other.pVertex_start and self.pVertex_end == other.pVertex_end

    def tojson(self) -> str:
        """Convert an edge object to a json string."""
        return json.dumps(self, cls=EdgeClassEncoder, sort_keys=True, indent=4)