import json
from json import JSONEncoder
import numpy as np
from typing import List
from pyearth.toolbox.mesh.vertex import pyvertex

class CircleClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pyvertex):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)

class pycircle(object):
    """The pycircle class represents a circle on a sphere."""

    def __init__(self, pVertex_center_in: pyvertex, aVertex_circle_in: List[pyvertex]):
        """
        Initialize a pycircle object.

        Args:
            pVertex_center_in (pyvertex): The center vertex of the circle.
            aVertex_circle_in (List[pyvertex]): A list of vertices that form the circle.
        """
        self.pVertex_center: pyvertex = pVertex_center_in
        self.aVertex_circle: List[pyvertex] = aVertex_circle_in

    def tojson(self) -> str:
        """
        Convert the pycircle object to a JSON string.

        Returns:
            str: A JSON string representing the pycircle object.
        """
        return json.dumps(self, cls=CircleClassEncoder, sort_keys=True, indent=4)
