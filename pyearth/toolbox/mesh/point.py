import json
from json import JSONEncoder
import importlib.util
import numpy as np

from pyearth.gis.gdal.write.vector.gdal_export_point_to_vector_file import (
    export_point_as_polygon_file,
)

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyearth.gis.geometry.kernel import (
        calculate_distance_based_on_longitude_latitude,
    )
else:
    from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import (
        calculate_distance_based_on_longitude_latitude,
    )

iPrecision_default = 12  # used for comparison


class PointClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)


class pypoint(object):
    """
    The point class

    Args:
        object (_type_): None

    Returns:
        pypoint: A point object
    """

    def __init__(self, aParameter):
        """
        Initialize a point object.

        Args:
            aParameter (dict): A dictionary containing point parameters.
                               Expected keys: 'dLongitude_degree', 'dLatitude_degree'.
                               Optional keys: 'x', 'y', 'z', 'dElevation'.
        """
        self.dX_meter = float(aParameter.get("x", -9999.0))
        self.dY_meter = float(aParameter.get("y", -9999.0))
        self.dZ_meter = float(aParameter.get("z", -9999.0))
        self.dElevation = float(aParameter.get("dElevation", 0.0))

        # dLongitude and dLatitude are always required
        if (
            "dLongitude_degree" not in aParameter
            or "dLatitude_degree" not in aParameter
        ):
            raise ValueError(
                "Initialization of pypoint failed: 'dLongitude_degree' and 'dLatitude_degree' are required."
            )

        self.dLongitude_degree = float(aParameter["dLongitude_degree"])
        self.dLatitude_degree = float(aParameter["dLatitude_degree"])

        self.dLongitude_radian = np.radians(self.dLongitude_degree)
        self.dLatitude_radian = np.radians(self.dLatitude_degree)

        # Recalculate x, y, z based on dLongitude and dLatitude if default values are used
        if (
            self.dX_meter == -9999.0
            and self.dY_meter == -9999.0
            and self.dZ_meter == -9999.0
        ):
            self.dX_meter, self.dY_meter, self.dZ_meter = self.calculate_xyz()

        self.wkt = self.towkt()

        return

    def __repr__(self):
        return f"pypoint(dLongitude_degree={self.dLongitude_degree}, dLatitude_degree={self.dLatitude_degree}, dElevation={self.dElevation})"

    def __str__(self):
        return f"({self.dLongitude_degree}, {self.dLatitude_degree}, {self.dElevation})"

    def toNvector(self):
        from pyearth.toolbox.mesh.nvector import pynvector

        """
        Note: replicated in LatLon_NvectorEllipsoidal

        Returns:
            pynvector: A nvector object
        """

        a = self.dLatitude_radian
        b = self.dLongitude_radian
        c = np.sin(a)
        e = np.cos(a)
        d = np.sin(b)
        f = np.cos(b)
        # // right-handed vector: x -> 0°E,0°N; y -> 90°E,0°N, z -> 90°N
        x = e * f
        y = e * d
        z = c
        point = dict()
        point["x"] = x
        point["y"] = y
        point["z"] = z
        pNvector = pynvector(point)
        return pNvector

    def __hash__(self):
        return hash(
            (
                round(self.dLongitude_degree, iPrecision_default),
                round(self.dLatitude_degree, iPrecision_default),
            )
        )

    def __eq__(self, other):
        """
        Check whether two points are equivalent.
        Uses numpy.isclose for robust floating-point comparison.
        """
        if not isinstance(other, pypoint):
            return NotImplemented

        dThreshold_in = 10.0 ** (-1 * iPrecision_default)
        return np.isclose(
            self.dLongitude_degree, other.dLongitude_degree, atol=dThreshold_in, rtol=0
        ) and np.isclose(
            self.dLatitude_degree, other.dLatitude_degree, atol=dThreshold_in, rtol=0
        )

    def __ne__(self, other):
        """
        Check whether two points are equivalent

        Args:
            other (pypoint): The other point

        Returns:
            int: 0 if equivalent, 1 if not
        """
        return not self.__eq__(other)

    def calculate_distance(self, other):
        """
        Calculate the distance between two points

        Args:
            other (pypoint): The other point

        Returns:
            float: The great circle distance
        """
        dDistance = 0.0
        lon1 = self.dLongitude_degree
        lat1 = self.dLatitude_degree
        lon2 = other.dLongitude_degree
        lat2 = other.dLatitude_degree
        dDistance = calculate_distance_based_on_longitude_latitude(
            lon1, lat1, lon2, lat2
        )
        return dDistance

    def calculate_buffer_zone_point(self, dRadius, dBearing=90):
        # Create a geodesic object
        from geographiclib.geodesic import Geodesic
        geod = Geodesic.WGS84  # the default is WGS84
        # Calculate the geodesic buffer
        pPoint_buffer = geod.Direct(
            self.dLatitude_degree, self.dLongitude_degree, dBearing, dRadius
        )
        # Extract the latitude and longitude of the buffer point
        # create a point object using the buffer point
        point0 = dict()
        point0["dLongitude_degree"] = pPoint_buffer["lon2"]
        point0["dLatitude_degree"] = pPoint_buffer["lat2"]
        pPoint_out = pypoint(point0)
        return pPoint_out

    def calculate_buffer_zone_circle(self, dRadius, nPoint=360, sFilename_out=None):
        # Create a geodesic object
        from geographiclib.geodesic import Geodesic
        geod = Geodesic.WGS84  # the default is WGS84
        aPoint = []
        # Calculate the geodesic buffer
        for i in range(0, 360, 360 // nPoint):
            pPoint_buffer = geod.Direct(
                self.dLatitude_degree, self.dLongitude_degree, i, dRadius
            )
            point0 = dict()
            point0["dLongitude_degree"] = pPoint_buffer["lon2"]
            point0["dLatitude_degree"] = pPoint_buffer["lat2"]
            pPoint_out = pypoint(point0)
            aPoint.append(pPoint_out)

        if sFilename_out is not None:
            # save as a geojson file
            export_point_as_polygon_file(aPoint, sFilename_out)

        return aPoint

    def calculate_xyz(self):
        """
        Calculate the x, y, z based on dLongitude and dLatitude

        Returns:
            tuple: The x, y, z
        """
        dX_meter = 0.0
        dY_meter = 0.0
        dZ_meter = 0.0
        dRadius = 6371000.0  # earth radius in meter
        dX_meter = (
            dRadius * np.cos(self.dLatitude_radian) * np.cos(self.dLongitude_radian)
        )
        dY_meter = (
            dRadius * np.cos(self.dLatitude_radian) * np.sin(self.dLongitude_radian)
        )
        dZ_meter = dRadius * np.sin(self.dLatitude_radian)
        return dX_meter, dY_meter, dZ_meter

    def tojson(self):
        """
        Convert a point object to a json string

        Returns:
            json str: A json string
        """
        aSkip = ["dLongitude_radian", "dLatitude_radian"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        # sJson = json.dumps(self.__dict__, \
        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=PointClassEncoder
        )
        return sJson

    def towkt(self):
        """
        Convert a point object to a WKT string

        Returns:
            str: A WKT string
        """
        sWKT = "POINT ("
        sWKT += str(self.dLongitude_degree) + " "
        sWKT += str(self.dLatitude_degree) + ")"
        self.wkt = sWKT
        return sWKT
