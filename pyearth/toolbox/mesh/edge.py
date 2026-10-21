"""
PyEdge class for representing directed edges connecting vertices on a sphere.

Extends the pyline class from pyearth to add topological edge attributes
(lEdgeID, lEdgeIndex, lIndex_upstream, lIndex_downstream) and functionality.
"""

import json
from json import JSONEncoder
from typing import Any, Dict, Optional, Tuple

import numpy as np

from pyearth.toolbox.mesh.line import pyline
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.vertex import pyvertex


class EdgeClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyedge objects.

    Handles numpy data types and pyvertex objects, converting them to
    native Python types for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            return obj
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)


class pyedge(pyline):
    """
    Edge class for topological network and mesh representations.

    Extends the pyline class to include edge-specific topological attributes.
    An edge represents a directed connection between two pyvertex nodes.

    Attributes:
        lEdgeID (int): Unique identifier for the edge (default: -1)
        lEdgeIndex (int): Array index for the edge (default: -1)
        lIndex_upstream (int): Index of upstream connected edge (default: -1)
        lIndex_downstream (int): Index of downstream connected edge (default: -1)
        pVertex_start (pyvertex): Starting vertex of the edge
        pVertex_end (pyvertex): Ending vertex of the edge

    Inherits from pyline:
        pPoint_start (pypoint): Start point of the line
        pPoint_end (pypoint): End point of the line
        dLength (float): Length of the edge in meters
        dSlope (float): Slope of the edge
    """

    def __init__(self, pVertex_start_in, pVertex_end_in) -> None:
        """
        Initialize a pyedge object.

        Args:
            pVertex_start_in (pyvertex or pypoint): The starting vertex
            pVertex_end_in (pyvertex or pypoint): The ending vertex

        Raises:
            TypeError: If inputs are not pyvertex or pypoint objects
            ValueError: If start and end vertices are identical (zero length)
        """
        if pVertex_start_in is None or pVertex_end_in is None:
            raise ValueError("Start and end vertices cannot be None")

        if hasattr(pVertex_start_in, "__dict__") and hasattr(
            pVertex_end_in, "__dict__"
        ):
            pPoint_start_in = pypoint(pVertex_start_in.__dict__)
            pPoint_end_in = pypoint(pVertex_end_in.__dict__)
        else:
            raise TypeError(
                "Start and end vertices must be pyvertex or pypoint objects"
            )

        if pPoint_start_in == pPoint_end_in:
            raise ValueError(
                "Start and end vertices cannot be identical (zero-length edge)"
            )

        super().__init__(pPoint_start_in, pPoint_end_in)

        # Edge-specific attributes with defaults
        self.lEdgeID: int = -1
        self.lEdgeIndex: int = -1
        self.lIndex_upstream: int = -1
        self.lIndex_downstream: int = -1

        self.pVertex_start = pVertex_start_in
        self.pVertex_end = pVertex_end_in

    def __repr__(self) -> str:
        return (
            f"pyedge(ID={self.lEdgeID}, Index={self.lEdgeIndex}, "
            f"Length={self.dLength:.2f}m, "
            f"Upstream={self.lIndex_upstream}, Downstream={self.lIndex_downstream})"
        )

    def __str__(self) -> str:
        return f"pyedge(ID={self.lEdgeID}, Length={self.dLength:.2f}m)"

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, pyedge):
            return NotImplemented
        return super().__eq__(other)

    def tojson(self) -> str:
        """Convert edge object to a JSON string."""
        aSkip = [
            "dLongitude_radian",
            "dLatitude_radian",
            "wkt",
            "pPoint_start",
            "pPoint_end",
        ]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=EdgeClassEncoder
        )
        return sJson

    def set_edge_id(self, lEdgeID: int) -> None:
        """Set the edge ID."""
        if not isinstance(lEdgeID, (int, np.integer)):
            raise TypeError(f"Edge ID must be an integer, got {type(lEdgeID)}")
        self.lEdgeID = int(lEdgeID)

    def set_edge_index(self, lEdgeIndex: int) -> None:
        """Set the edge array index."""
        if not isinstance(lEdgeIndex, (int, np.integer)):
            raise TypeError(f"Edge index must be an integer, got {type(lEdgeIndex)}")
        self.lEdgeIndex = int(lEdgeIndex)

    def set_upstream_index(self, lIndex_upstream: int) -> None:
        """Set the upstream edge index."""
        if not isinstance(lIndex_upstream, (int, np.integer)):
            raise TypeError(
                f"Upstream index must be an integer, got {type(lIndex_upstream)}"
            )
        self.lIndex_upstream = int(lIndex_upstream)

    def set_downstream_index(self, lIndex_downstream: int) -> None:
        """Set the downstream edge index."""
        if not isinstance(lIndex_downstream, (int, np.integer)):
            raise TypeError(
                f"Downstream index must be an integer, got {type(lIndex_downstream)}"
            )
        self.lIndex_downstream = int(lIndex_downstream)

    def is_valid(self) -> bool:
        """Check if the edge has valid attributes."""
        has_valid_length = self.dLength > 0
        has_valid_vertices = (
            hasattr(self, "pVertex_start")
            and hasattr(self, "pVertex_end")
            and self.pVertex_start.is_valid()
            and self.pVertex_end.is_valid()
        )
        has_valid_id = self.lEdgeID > 0 or self.lEdgeIndex >= 0
        return has_valid_length and has_valid_vertices and has_valid_id

    def copy(self) -> "pyedge":
        """Create a deep copy of the edge."""
        new_edge = pyedge(self.pVertex_start.copy(), self.pVertex_end.copy())
        new_edge.lEdgeID = self.lEdgeID
        new_edge.lEdgeIndex = self.lEdgeIndex
        new_edge.lIndex_upstream = self.lIndex_upstream
        new_edge.lIndex_downstream = self.lIndex_downstream
        return new_edge

    def reverse(self) -> None:
        """
        Reverse the direction of the edge in place by swapping start and end vertices,
        and swapping upstream/downstream indices.
        """
        super().reverse()
        self.pVertex_start, self.pVertex_end = self.pVertex_end, self.pVertex_start
        self.lIndex_upstream, self.lIndex_downstream = (
            self.lIndex_downstream,
            self.lIndex_upstream,
        )

    def reversed(self) -> "pyedge":
        """
        Create a reversed copy of the edge with swapped start and end vertices,
        and swapped upstream/downstream indices.
        """
        reversed_edge = self.copy()
        reversed_edge.reverse()
        return reversed_edge

    def get_midpoint(self) -> pyvertex:
        """Calculate and return the midpoint of the edge as a pyvertex."""
        mid_lon = (
            self.pVertex_start.dLongitude_degree + self.pVertex_end.dLongitude_degree
        ) / 2.0
        mid_lat = (
            self.pVertex_start.dLatitude_degree + self.pVertex_end.dLatitude_degree
        ) / 2.0

        param = {"dLongitude_degree": mid_lon, "dLatitude_degree": mid_lat}
        return pyvertex(param)

    def is_connected_to(self, other: "pyedge") -> bool:
        """Check if this edge shares a vertex with another edge."""
        if not isinstance(other, pyedge):
            return False

        return (
            self.pVertex_end == other.pVertex_start
            or self.pVertex_end == other.pVertex_end
            or self.pVertex_start == other.pVertex_start
            or self.pVertex_start == other.pVertex_end
        )

    def is_upstream_of(self, other: "pyedge") -> bool:
        """Check if this edge is directly upstream of another edge."""
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_end == other.pVertex_start

    def is_downstream_of(self, other: "pyedge") -> bool:
        """Check if this edge is directly downstream of another edge."""
        if not isinstance(other, pyedge):
            return False
        return self.pVertex_start == other.pVertex_end

    def check_vertex_on_edge(self, pVertex_in: pyvertex) -> Tuple[int, float, float]:
        """
        Check if a vertex lies on this edge.

        Calculates whether the sum of distances from the vertex to both endpoints
        equals the edge length within a tolerance.

        Args:
            pVertex_in (pyvertex): The vertex to check.

        Returns:
            Tuple[int, float, float]:
                - int: 1 if vertex is on edge, 0 otherwise
                - float: Distance from start vertex to the test vertex (meters)
                - float: Difference between (d_start + d_end) and edge length (meters)
        """
        d_start = self.pVertex_start.calculate_distance(pVertex_in)
        d_end = pVertex_in.calculate_distance(self.pVertex_end)
        diff = abs((d_start + d_end) - self.dLength)

        # Tolerant threshold: within 1 meter or 0.01% of edge length
        if diff < 1.0 or (self.dLength > 0 and diff / self.dLength < 1e-4):
            return 1, d_start, diff
        return 0, d_start, diff
