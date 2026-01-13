import json
from json import JSONEncoder
import numpy as np
from typing import List
from pyearth.toolbox.mesh.point import pypoint


class CircleClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pypoint):
            return json.loads(obj.tojson())
        return JSONEncoder.default(self, obj)


class pycircle(object):
    """The pycircle class represents a circle on a sphere."""

    def __init__(self, pPoint_center_in: pypoint, aPoint_circle_in: List[pypoint]):
        """
        Initialize a pycircle object.

        Args:
            pPoint_center_in (pypoint): The center point of the circle.
            aPoint_circle_in (List[pypoint]): A list of points that form the circle.
        """
        self.pPoint_center: pypoint = pPoint_center_in
        self.aPoint_circle: List[pypoint] = aPoint_circle_in

    def tojson(self) -> str:
        """
        Convert the pycircle object to a JSON string.

        Returns:
            str: A JSON string representing the pycircle object.
        """
        return json.dumps(self, cls=CircleClassEncoder, sort_keys=True, indent=4)
