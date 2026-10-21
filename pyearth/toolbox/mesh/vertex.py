"""
PyVertex class for representing topological vertices/points on a sphere.

Extends the pypoint class from pyearth to add topological vertex attributes
(lVertexID, lVertexIndex, lFlowlineID) and functionality.
"""

import json
from json import JSONEncoder
from typing import Any, Dict, Optional

import numpy as np

from pyearth.toolbox.mesh.point import pypoint


class VertexClassEncoder(JSONEncoder):
    """
    Custom JSON encoder for pyvertex objects.

    Handles numpy data types and converts them to native Python types
    for JSON serialization.
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)


class pyvertex(pypoint):
    """
    Vertex class for topological network and mesh representations.

    Extends the pypoint class to include vertex-specific topological attributes.
    A vertex represents a node that may be shared between multiple edges/lines,
    flowlines, or mesh cells.

    Attributes:
        lVertexIndex (int): Array index for the vertex (default: -1)
        lVertexID (int): Unique identifier for the vertex (default: 1)
        lFlowlineID (int): Associated flowline or feature ID (default: -1)

    Inherits from pypoint:
        dLongitude_degree (float): Longitude in degrees
        dLatitude_degree (float): Latitude in degrees
        dX_meter (float): X coordinate in meters (ECEF)
        dY_meter (float): Y coordinate in meters (ECEF)
        dZ_meter (float): Z coordinate in meters (ECEF)
        dElevation (float): Elevation value
    """

    def __init__(self, aParameter: Dict[str, Any]) -> None:
        """
        Initialize a vertex object.

        Args:
            aParameter (dict): Dictionary containing vertex parameters.
                             Must include 'dLongitude_degree' and 'dLatitude_degree'.

        Raises:
            ValueError: If required parameters are missing
            TypeError: If aParameter is not a dictionary
        """
        if not isinstance(aParameter, dict):
            raise TypeError(f"aParameter must be a dictionary, got {type(aParameter)}")

        # Initialize parent class (pypoint)
        super().__init__(aParameter)

        # Vertex-specific topological attributes with defaults
        self.lVertexIndex: int = int(aParameter.get("lVertexIndex", -1))
        self.lVertexID: int = int(aParameter.get("lVertexID", 1))
        self.lFlowlineID: int = int(aParameter.get("lFlowlineID", -1))

    def __repr__(self) -> str:
        return (
            f"pyvertex(ID={self.lVertexID}, Index={self.lVertexIndex}, "
            f"Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f}, "
            f"FlowlineID={self.lFlowlineID})"
        )

    def __str__(self) -> str:
        return f"pyvertex(ID={self.lVertexID}, Lon={self.dLongitude_degree:.6f}, Lat={self.dLatitude_degree:.6f})"

    def __hash__(self) -> int:
        return super().__hash__()

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, pyvertex):
            return NotImplemented
        return super().__eq__(other)

    def tojson(self) -> str:
        """
        Convert vertex object to a JSON string.

        Serializes all vertex attributes except internal computed values.
        """
        aSkip = ["dLongitude_radian", "dLatitude_radian", "wkt"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=VertexClassEncoder
        )
        return sJson

    def set_vertex_id(self, lVertexID: int) -> None:
        """Set the vertex ID."""
        if not isinstance(lVertexID, (int, np.integer)):
            raise TypeError(f"Vertex ID must be an integer, got {type(lVertexID)}")
        self.lVertexID = int(lVertexID)

    def set_vertex_index(self, lVertexIndex: int) -> None:
        """Set the vertex array index."""
        if not isinstance(lVertexIndex, (int, np.integer)):
            raise TypeError(
                f"Vertex index must be an integer, got {type(lVertexIndex)}"
            )
        self.lVertexIndex = int(lVertexIndex)

    def set_flowline_id(self, lFlowlineID: int) -> None:
        """Set the associated flowline or feature ID."""
        if not isinstance(lFlowlineID, (int, np.integer)):
            raise TypeError(f"Flowline ID must be an integer, got {type(lFlowlineID)}")
        self.lFlowlineID = int(lFlowlineID)

    def toNvector(self, use_high_precision: bool = False):
        """
        Convert geographic coordinates to n-vector representation.

        Args:
            use_high_precision (bool): Use high precision float for calculations.

        Returns:
            pynvector: An n-vector object with unit length
        """
        from pyearth.toolbox.mesh.nvector import pynvector

        if use_high_precision:
            dtype = getattr(np, "float128", np.longdouble)
            a = dtype(self.dLatitude_radian)
            b = dtype(self.dLongitude_radian)
        else:
            a = self.dLatitude_radian
            b = self.dLongitude_radian

        c = np.sin(a)
        e = np.cos(a)
        d = np.sin(b)
        f = np.cos(b)

        point = {
            "x": e * f,
            "y": e * d,
            "z": c,
        }
        return pynvector(point, use_high_precision=use_high_precision)

    def is_valid(self) -> bool:
        """Check if the vertex has valid attributes."""
        has_valid_coords = (
            -180 <= self.dLongitude_degree <= 180 and -90 <= self.dLatitude_degree <= 90
        )
        has_valid_id = self.lVertexID > 0 or self.lVertexIndex >= 0
        return has_valid_coords and has_valid_id

    def copy(self) -> "pyvertex":
        """Create a deep copy of the vertex."""
        param = {
            "dLongitude_degree": self.dLongitude_degree,
            "dLatitude_degree": self.dLatitude_degree,
            "x": self.dX_meter,
            "y": self.dY_meter,
            "z": self.dZ_meter,
            "dElevation": self.dElevation,
            "lVertexIndex": self.lVertexIndex,
            "lVertexID": self.lVertexID,
            "lFlowlineID": self.lFlowlineID,
        }
        return pyvertex(param)
