import json
from json import JSONEncoder
from typing import List, Tuple, Optional
import numpy as np
from osgeo import ogr
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline
from pyearth.toolbox.mesh.polyline import pypolyline
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area


class PolygonClassEncoder(JSONEncoder):
    """Custom JSON encoder for pypolygon objects and their dependencies."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, list):
            pass
        if isinstance(obj, pypoint):
            return json.loads(obj.tojson())
        if isinstance(obj, pyline):
            return obj.lLineID
        if isinstance(obj, pypolyline):
            return obj.lLineID
        if isinstance(obj, pypolygon):
            return obj.lPolygonID
        return super().default(obj)


class pypolygon:
    """
    Represents a polygon composed of connected line segments forming a closed shape.

    A polygon is a closed geometric figure made up of line segments where the first
    and last points are identical. This class provides geometric operations including
    area calculations, boundary checks, and topological relationships between polygons.

    Attributes:
        aLine (List[pyline]): List of line segments forming the polygon boundary.
        aPoint (List[pypoint]): List of vertices (first and last are identical).
        pPoint_center (pypoint): Centroid of the polygon.
        dLongitude_center_degree (float): Longitude of the center in degrees.
        dLatitude_center_degree (float): Latitude of the center in degrees.
        nLine (int): Number of line segments.
        nPoint (int): Number of vertices.
        dLength (float): Characteristic length (sqrt of area).
        dArea (float): Area of the polygon in square meters.
        pBound (Tuple[float, float, float, float]): Bounding box (lon_min, lat_min, lon_max, lat_max).
        wkt (str): WKT representation of the polygon.
        lPolygonID (int): Optional identifier (default: -1).
        dX_center_meter (float): X coordinate of center (default: 0.0).
        dY_center_meter (float): Y coordinate of center (default: 0.0).
        dz_center (float): Z coordinate of center (default: 0.0).

    Example:
        >>> from pyearth.toolbox.mesh.point import pypoint
        >>> from pyearth.toolbox.mesh.line import pyline
        >>> p1 = pypoint({'dLongitude_degree': -122.4, 'dLatitude_degree': 37.8})
        >>> p2 = pypoint({'dLongitude_degree': -122.3, 'dLatitude_degree': 37.8})
        >>> p3 = pypoint({'dLongitude_degree': -122.3, 'dLatitude_degree': 37.7})
        >>> p4 = pypoint({'dLongitude_degree': -122.4, 'dLatitude_degree': 37.7})
        >>> lines = [pyline(p1, p2), pyline(p2, p3), pyline(p3, p4), pyline(p4, p1)]
        >>> points = [p1, p2, p3, p4, p1]
        >>> polygon = pypolygon(-122.35, 37.75, lines, points)
        >>> print(f"Area: {polygon.calculate_polygon_area():.2f} m²")
    """

    def __init__(
        self, dLon: float, dLat: float, aLine: List[pyline], aPoint: List[pypoint]
    ):
        """
        Initialize a polygon object.

        Args:
            dLon (float): The longitude of the center in degrees.
            dLat (float): The latitude of the center in degrees.
            aLine (List[pyline]): A list of line segments that define the polygon boundary.
            aPoint (List[pypoint]): A list of vertices that define the polygon (first and last must be identical).

        Raises:
            ValueError: If polygon has fewer than 3 lines, or if first and last points don't match.
            TypeError: If inputs are not of expected types.
        """
        if not isinstance(aLine, list) or not isinstance(aPoint, list):
            raise TypeError("aLine and aPoint must be lists.")

        nLine = len(aLine)
        if nLine < 3:
            raise ValueError("Polygon must have at least 3 line segments.")

        if not all(isinstance(line, pyline) for line in aLine):
            raise TypeError("All elements in aLine must be pyline objects.")

        if not all(isinstance(point, pypoint) for point in aPoint):
            raise TypeError("All elements in aPoint must be pypoint objects.")

        # Verify polygon is closed (first and last points are identical)
        if len(aPoint) > 0 and aPoint[0] != aPoint[-1]:
            raise ValueError(
                "Polygon must be closed: first and last points must be identical."
            )

        # Initialize attributes
        self.lPolygonID: int = -1
        self.nPoint: int = len(aPoint)
        self.nLine: int = nLine
        self.dLength: float = 0.0
        self.dArea: float = 0.0
        self.dX_center_meter: float = 0.0
        self.dY_center_meter: float = 0.0
        self.dz_center: float = 0.0
        self.dLongitude_center_degree: float = float(dLon)
        self.dLatitude_center_degree: float = float(dLat)

        self.aLine: List[pyline] = aLine
        self.aPoint: List[pypoint] = aPoint

        # Create center point
        pPoint = {
            "dLongitude_degree": self.dLongitude_center_degree,
            "dLatitude_degree": self.dLatitude_center_degree,
        }
        self.pPoint_center: pypoint = pypoint(pPoint)

        # Calculate derived properties
        self.pBound: Tuple[float, float, float, float] = self.calculate_polygon_bound()
        self.wkt: str = self.towkt()

    def calculate_polygon_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the polygon.

        Returns:
            Tuple[float, float, float, float]: Bounding box (lon_min, lat_min, lon_max, lat_max) in degrees.
        """
        if self.nPoint == 0:
            return (0.0, 0.0, 0.0, 0.0)

        # Use numpy for efficient calculation
        lons = np.array([p.dLongitude_degree for p in self.aPoint])
        lats = np.array([p.dLatitude_degree for p in self.aPoint])

        dLon_min = float(np.min(lons))
        dLon_max = float(np.max(lons))
        dLat_min = float(np.min(lats))
        dLat_max = float(np.max(lats))

        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def has_this_line(self, pLine_in: pyline) -> bool:
        """
        Check whether the polygon contains a specific line.

        Args:
            pLine_in (pyline): The line to be checked.

        Returns:
            bool: True if the line is found, False otherwise.
        """
        return any(pLine.is_overlap(pLine_in) for pLine in self.aLine)

    def which_line_cross_this_point(
        self, pPoint_in: pypoint
    ) -> Tuple[bool, Optional[pyline]]:
        """
        Find which line (edge) of the polygon contains a given point.

        Args:
            pPoint_in (pypoint): The point to be checked.

        Returns:
            Tuple[bool, Optional[pyline]]:
                - True if found with the line object
                - False if not found, with None
        """
        for pLine in self.aLine:
            iFlag, dummy, diff = pLine.check_point_on_line(pPoint_in)
            if iFlag:
                return True, pLine

        return False, None

    def calculate_polygon_area(self) -> float:
        """
        Calculate the area of the polygon using spherical geometry.

        Returns:
            float: The area in square meters.

        Raises:
            ValueError: If polygon has fewer than 3 points.
        """
        if self.nPoint < 4:  # Need at least 3 unique points + closing point
            raise ValueError(
                "Polygon must have at least 3 unique vertices to calculate area."
            )

        lons = [p.dLongitude_degree for p in self.aPoint]
        lats = [p.dLatitude_degree for p in self.aPoint]

        self.dArea = calculate_polygon_area(lons, lats)
        return self.dArea

    def calculate_line_length(self) -> float:
        """
        Calculate the characteristic length of the polygon.

        The characteristic length is defined as the square root of the area,
        providing a single length metric for the polygon.

        Returns:
            float: The characteristic length in meters.
        """
        if self.dArea == 0.0:
            self.calculate_polygon_area()

        dLength_line = np.sqrt(self.dArea)
        self.dLength = dLength_line
        return dLength_line

    def share_line(self, other: "pypolygon") -> bool:
        """
        Check whether this polygon shares a line (edge) with another polygon.

        Args:
            other (pypolygon): The other polygon to compare with.

        Returns:
            bool: True if polygons share at least one edge, False otherwise.

        Raises:
            TypeError: If other is not a pypolygon object.
        """
        if not isinstance(other, pypolygon):
            raise TypeError("Argument must be a pypolygon object.")

        for pLine in self.aLine:
            for pLine2 in other.aLine:
                if pLine.is_overlap(pLine2):
                    return True

        return False

    def __eq__(self, other: object) -> bool:
        """
        Check if two polygons are equivalent.

        Two polygons are considered equal if they have the same center coordinates.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if polygons have the same center, False otherwise.
        """
        if not isinstance(other, pypolygon):
            return NotImplemented

        return self.pPoint_center == other.pPoint_center and len(self.aLine) == len(
            other.aLine
        )

    def __ne__(self, other: object) -> bool:
        """
        Check if two polygons are not equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if polygons are not equivalent.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self) -> int:
        """
        Generate hash for the polygon based on its center point.

        Returns:
            int: Hash value for the polygon.
        """
        return hash((self.pPoint_center, self.nLine))

    def __repr__(self) -> str:
        """
        Return a detailed string representation for debugging.

        Returns:
            str: Developer-friendly representation.
        """
        return (
            f"pypolygon(nLine={self.nLine}, nPoint={self.nPoint}, "
            f"area={self.dArea:.2f}m², "
            f"center=({self.dLongitude_center_degree:.6f}°, "
            f"{self.dLatitude_center_degree:.6f}°), "
            f"ID={self.lPolygonID})"
        )

    def __str__(self) -> str:
        """
        Return a user-friendly string representation.

        Returns:
            str: Human-readable description.
        """
        return (
            f"Polygon with {self.nLine} edges, {self.nPoint} vertices at "
            f"({self.dLongitude_center_degree:.4f}°, "
            f"{self.dLatitude_center_degree:.4f}°)"
        )

    def tojson(self) -> str:
        """
        Convert the polygon object to a JSON string.

        Returns:
            str: JSON representation of the polygon.

        Note:
            The 'aLine' attribute is excluded from JSON output to avoid
            circular references and reduce size.
        """
        aSkip = ["aLine"]
        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=PolygonClassEncoder
        )
        return sJson

    def towkt(self) -> str:
        """
        Convert the polygon object to a WKT (Well-Known Text) string.

        Returns:
            str: WKT representation of the polygon.

        Raises:
            ValueError: If polygon has fewer than 3 unique points.
        """
        if self.nPoint < 4:
            raise ValueError(
                "Polygon must have at least 3 unique vertices for WKT conversion."
            )

        pGeometry = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)

        for point in self.aPoint:
            ring.AddPoint(point.dLongitude_degree, point.dLatitude_degree)

        ring.CloseRings()
        pGeometry.AddGeometry(ring)

        sWKT = pGeometry.ExportToWkt()
        self.wkt = sWKT
        return sWKT

    def update_wkt(self) -> str:
        """
        Update and return the WKT string of the polygon.

        Returns:
            str: Updated WKT representation.
        """
        return self.towkt()
