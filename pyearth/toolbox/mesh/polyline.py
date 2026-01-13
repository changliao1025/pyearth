import os
import copy
import json
import importlib.util
from json import JSONEncoder
from typing import List, Tuple, Optional
import numpy as np
from osgeo import ogr
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.line import pyline
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


class PolylineClassEncoder(JSONEncoder):
    """Custom JSON encoder for pypolyline objects and their dependencies."""

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

        return JSONEncoder.default(self, obj)


class pypolyline:
    """
    Represents a polyline composed of connected line segments.

    A polyline is a continuous path made up of one or more line segments,
    where each segment connects to the next. This class provides geometric
    operations including length calculations, splitting, reversing, and
    buffer zone generation.

    Attributes:
        aLine (List[pyline]): List of line segments composing the polyline.
        aPoint (List[pypoint]): List of all points along the polyline.
        pPoint_start (pypoint): Starting point of the polyline.
        pPoint_end (pypoint): Ending point of the polyline.
        nLine (int): Number of line segments.
        nPoint (int): Number of points (nLine + 1).
        dLength (float): Total geodesic length in meters.
        wkt (str): WKT representation of the polyline.
        pBound (Tuple[float, float, float, float]): Bounding box.
        lLineID (int): Optional identifier (default: -1).
        lLineIndex (int): Optional index (default: -1).
        iFlag_right (int): Right flag (default: 0).
        iFlag_left (int): Left flag (default: 0).
        dSinuosity (float): Sinuosity value (set by calculate_polyline_sinuosity).
    """

    def __init__(self, aLine: List[pyline]):
        """
        Initialize a polyline object.

        Args:
            aLine (List[pyline]): A list of connected line segments.

        Raises:
            ValueError: If aLine is empty or lines are not connected.
            TypeError: If aLine contains non-pyline objects.
        """
        if not aLine:
            raise ValueError("Cannot create polyline from empty line list.")

        if not all(isinstance(line, pyline) for line in aLine):
            raise TypeError("All elements in aLine must be pyline objects.")

        # Validate that lines are connected
        for i in range(len(aLine) - 1):
            if aLine[i].pPoint_end != aLine[i + 1].pPoint_start:
                raise ValueError(f"Lines are not connected at index {i}.")

        self.lLineID: int = -1
        self.lLineIndex: int = -1
        self.iFlag_right: int = 0
        self.iFlag_left: int = 0

        # Store the line segments (CRITICAL FIX)
        self.aLine: List[pyline] = aLine
        self.nLine: int = len(aLine)

        # Set start and end points
        self.pPoint_start: pypoint = aLine[0].pPoint_start
        self.pPoint_end: pypoint = aLine[self.nLine - 1].pPoint_end

        # Build point list from line segments
        self.aPoint: List[pypoint] = []
        for i in range(self.nLine):
            self.aPoint.append(aLine[i].pPoint_start)
        self.aPoint.append(aLine[self.nLine - 1].pPoint_end)
        self.nPoint: int = len(self.aPoint)

        # Calculate derived properties
        self.dLength: float = self.calculate_length()
        self.wkt: str = self.towkt()
        self.pBound: Tuple[float, float, float, float] = self.calculate_line_bound()

    def __hash__(self) -> int:
        """
        Generate hash for the polyline based on all its points.

        Returns:
            int: Hash value for the polyline.
        """
        return hash(tuple(self.aPoint))

    def __repr__(self) -> str:
        """
        Return a detailed string representation for debugging.

        Returns:
            str: Developer-friendly representation.
        """
        return (
            f"pypolyline(nLine={self.nLine}, nPoint={self.nPoint}, "
            f"length={self.dLength:.2f}m, "
            f"start=({self.pPoint_start.dLongitude_degree:.6f}, "
            f"{self.pPoint_start.dLatitude_degree:.6f}), "
            f"end=({self.pPoint_end.dLongitude_degree:.6f}, "
            f"{self.pPoint_end.dLatitude_degree:.6f}))"
        )

    def __str__(self) -> str:
        """
        Return a user-friendly string representation.

        Returns:
            str: Human-readable description.
        """
        return (
            f"Polyline with {self.nLine} segments, {self.nPoint} points, "
            f"length {self.dLength:.2f}m from "
            f"({self.pPoint_start.dLongitude_degree:.4f}°, "
            f"{self.pPoint_start.dLatitude_degree:.4f}°) to "
            f"({self.pPoint_end.dLongitude_degree:.4f}°, "
            f"{self.pPoint_end.dLatitude_degree:.4f}°)"
        )

    def calculate_length(self) -> float:
        """
        Calculate the total length of the polyline.

        Returns:
            float: The total geodesic length in meters.
        """
        self.dLength = sum(line.dLength for line in self.aLine)
        return self.dLength

    def calculate_line_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the polyline.

        Returns:
            Tuple[float, float, float, float]: Bounding box (lon_min, lat_min, lon_max, lat_max).

        Raises:
            ValueError: If WKT geometry is invalid.
        """
        pGeometry = ogr.CreateGeometryFromWkt(self.wkt)
        if pGeometry is None:
            raise ValueError("Invalid geometry in polyline")
        (minX, maxX, minY, maxY) = pGeometry.GetEnvelope()
        self.pBound = (float(minX), float(minY), float(maxX), float(maxY))
        return self.pBound

    def reverse(self) -> None:
        """
        Reverse the direction of the polyline.

        This method reverses both the point list and reconstructs the line
        segments in reverse order, updating start and end points accordingly.
        """
        # Reverse point list
        self.aPoint = list(reversed(self.aPoint))

        # Rebuild line segments from reversed points
        self.aLine = [
            pyline(self.aPoint[i], self.aPoint[i + 1])
            for i in range(len(self.aPoint) - 1)
        ]

        # Update endpoints
        self.pPoint_start = self.aLine[0].pPoint_start
        self.pPoint_end = self.aLine[-1].pPoint_end

    def split_line_by_length(self, dDistance: float) -> "pypolyline":
        """
        Split line segments that exceed the length threshold.

        Args:
            dDistance (float): Maximum length for each segment in meters.

        Returns:
            pypolyline: New polyline with subdivided segments.

        Raises:
            ValueError: If dDistance is not positive.
        """
        if dDistance <= 0:
            raise ValueError("Distance threshold must be positive.")

        aLine = []
        for line in self.aLine:
            if line.dLength > dDistance:
                aLine.extend(line.split_by_length(dDistance))
            else:
                aLine.append(line)

        pPolyline_out = pypolyline(aLine)
        pPolyline_out.copy_attributes(self)
        return pPolyline_out

    def split_by_length(self, dDistance: float) -> List["pypolyline"]:
        """
        Split polyline into multiple polylines based on length threshold.

        Uses recursive bisection to split the polyline until all segments
        are shorter than the threshold.

        Args:
            dDistance (float): Maximum length for each polyline in meters.

        Returns:
            List[pypolyline]: List of polylines, each shorter than threshold.

        Raises:
            ValueError: If dDistance is not positive.
        """
        if dDistance <= 0:
            raise ValueError("Distance threshold must be positive.")

        if self.dLength <= dDistance:
            return [self]

        # Split at midpoint
        mid_idx = len(self.aLine) // 2
        aLine0 = self.aLine[:mid_idx]
        aLine1 = self.aLine[mid_idx:]

        pPolyline0 = pypolyline(aLine0)
        pPolyline1 = pypolyline(aLine1)
        pPolyline0.copy_attributes(self)
        pPolyline1.copy_attributes(self)

        # Recursively split if still too long
        aPolyline = []
        if pPolyline0.dLength > dDistance:
            aPolyline.extend(pPolyline0.split_by_length(dDistance))
        else:
            aPolyline.append(pPolyline0)

        if pPolyline1.dLength > dDistance:
            aPolyline.extend(pPolyline1.split_by_length(dDistance))
        else:
            aPolyline.append(pPolyline1)

        return aPolyline

    def copy_attributes(self, other: "pypolyline") -> None:
        """
        Copy attributes from another polyline.

        Args:
            other (pypolyline): Source polyline to copy attributes from.
        """
        self.iFlag_right = other.iFlag_right
        self.iFlag_left = other.iFlag_left

    def calculate_polyline_sinuosity(self) -> float:
        """
        Calculate the sinuosity of the polyline.

        Sinuosity is the ratio of the actual path length to the straight-line
        distance between endpoints. A value of 1.0 indicates a straight line,
        higher values indicate more meandering.

        Returns:
            float: The sinuosity ratio (>= 1.0).
        """
        dDistance_straight = self.pPoint_start.calculate_distance(self.pPoint_end)
        if dDistance_straight == 0:
            return 1.0
        self.dSinuosity = self.dLength / dDistance_straight
        return self.dSinuosity

    def calculate_distance_to_point(self, pPoint: pypoint) -> Tuple[float, pypoint]:
        """
        Calculate the minimum distance from a point to this polyline.

        This method checks distances to all vertices and line segments,
        returning the minimum distance and the closest point on the polyline.

        Args:
            pPoint (pypoint): The point to measure distance from.

        Returns:
            Tuple[float, pypoint]:
                - Minimum distance in meters
                - Closest point on the polyline
        """
        # Check distances to all vertices
        dDistance_min_point = float("inf")
        pPoint_min_point = None

        for point in self.aPoint:
            dDistance = pPoint.calculate_distance(point)
            if dDistance < dDistance_min_point:
                dDistance_min_point = dDistance
                pPoint_min_point = point

        # Check distances to all line segments (FIXED: was using non-existent nEdge/aEdge)
        dDistance_min_edge = float("inf")
        pPoint_min_edge = None

        for line in self.aLine:
            dDistance, pPoint_out_edge = line.calculate_distance_to_point(pPoint)
            if dDistance < dDistance_min_edge:
                dDistance_min_edge = dDistance
                pPoint_min_edge = pPoint_out_edge

        # Return the minimum of vertex and edge distances
        if dDistance_min_point < dDistance_min_edge:
            return dDistance_min_point, pPoint_min_point
        else:
            return dDistance_min_edge, pPoint_min_edge

    def calculate_buffer_zone_polygon(
        self,
        dRadius: float,
        sFilename_out: Optional[str] = None,
        sFolder_out: Optional[str] = None,
    ) -> Tuple[str, List[pypoint], List]:
        """
        Calculate the buffer zone polygon around the polyline.

        Args:
            dRadius (float): Buffer radius in meters.
            sFilename_out (Optional[str]): Output file path for the buffer polygon.
            sFolder_out (Optional[str]): Output folder for intermediate buffer files.

        Returns:
            Tuple containing:
                - str: WKT representation of the buffer polygon
                - List[pypoint]: Vertices of the buffer polygon
                - List: Circle objects created

        Raises:
            ValueError: If dRadius is not positive.
        """
        if dRadius <= 0:
            raise ValueError("Buffer radius must be positive.")

        pMultiPolygon = ogr.Geometry(ogr.wkbMultiPolygon)
        aPoint_out = list()
        aCircle_out = list()
        for i in range(self.nLine):
            line = self.aLine[i]
            aPoint, aPoint_center, aPoint_circle, aCircle = (
                line.calculate_buffer_zone_polygon(dRadius)
            )
            aCircle_out.append(aCircle)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for pPoint in aPoint:
                ring.AddPoint(pPoint.dLongitude_degree, pPoint.dLatitude_degree)

            ring.CloseRings()
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            # Add the polygon to the MultiPolygon
            pMultiPolygon.AddGeometry(pPolygon)
            if sFolder_out is not None:
                sFilename_dummy = os.path.join(
                    sFolder_out, "buffer_zone_line_%d.geojson" % i
                )
                export_point_as_polygon_file(aPoint, sFilename_dummy)

        pUnionPolygon = pMultiPolygon.UnionCascaded()
        for i in range(pUnionPolygon.GetGeometryRef(0).GetPointCount()):
            lon, lat, _ = pUnionPolygon.GetGeometryRef(0).GetPoint(i)
            point2 = dict()
            point2["dLongitude_degree"] = lon
            point2["dLatitude_degree"] = lat
            pPoint2 = pypoint(point2)
            aPoint_out.append(pPoint2)

        if sFilename_out is not None:
            export_point_as_polygon_file(aPoint_out, sFilename_out)

        sWkt_buffer_polygon = pUnionPolygon.ExportToWkt()

        return sWkt_buffer_polygon, aPoint_out, aCircle_out

    def calculate_distance_to_polyline(self, pPolyline_other: "pypolyline") -> float:
        """
        Calculate the minimum distance between two polylines.

        Uses vectorized operations to efficiently compute all pairwise distances
        between points on both polylines.

        Args:
            pPolyline_other (pypolyline): The other polyline to measure distance to.

        Returns:
            float: Minimum distance in meters between the two polylines.
        """
        aPoint_a = np.array(
            [[v.dLongitude_degree, v.dLatitude_degree] for v in self.aPoint]
        )
        aPoint_b = np.array(
            [[v.dLongitude_degree, v.dLatitude_degree] for v in pPolyline_other.aPoint]
        )

        # Use broadcasting to calculate all pairwise distances at once
        lon1 = aPoint_a[:, 0:1]  # Convert to column vector
        lat1 = aPoint_a[:, 1:2]
        lon2 = aPoint_b[:, 0]  # Keep as row vector
        lat2 = aPoint_b[:, 1]
        # Vectorized distance calculation
        distances = calculate_distance_based_on_longitude_latitude(
            lon1, lat1, lon2, lat2
        )
        return np.min(distances)

    def calculate_bearing_angle(self) -> Optional[float]:
        """
        Calculate the bearing angle of the polyline.

        The bearing is calculated from the start point to the end point using
        the forward azimuth formula on a sphere.

        Returns:
            Optional[float]: Bearing angle in degrees (0-360°), where 0° is North,
                            or None if polyline has less than 2 points.
        """
        if self.nPoint < 2:
            return None

        dLon_start = self.pPoint_start.dLongitude_degree
        dLat_start = self.pPoint_start.dLatitude_degree
        dLon_end = self.pPoint_end.dLongitude_degree
        dLat_end = self.pPoint_end.dLatitude_degree

        # Convert to radians
        lat1_rad = np.radians(dLat_start)
        lat2_rad = np.radians(dLat_end)
        dlon_rad = np.radians(dLon_end - dLon_start)

        # Calculate bearing using the forward azimuth formula
        y = np.sin(dlon_rad) * np.cos(lat2_rad)
        x = np.cos(lat1_rad) * np.sin(lat2_rad) - np.sin(lat1_rad) * np.cos(
            lat2_rad
        ) * np.cos(dlon_rad)

        # Calculate bearing in radians
        bearing_rad = np.arctan2(y, x)

        # Convert to degrees and normalize to 0-360°
        bearing_deg = np.degrees(bearing_rad)
        bearing_deg = (bearing_deg + 360) % 360

        return bearing_deg

    def __eq__(self, other: object) -> bool:
        """
        Check whether two polylines are equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if polylines have identical line segments in order.
        """
        if not isinstance(other, pypolyline):
            return NotImplemented

        if len(self.aLine) != len(other.aLine):
            return False

        return all(line1 == line2 for line1, line2 in zip(self.aLine, other.aLine))

    def __ne__(self, other: object) -> bool:
        """
        Check whether two polylines are not equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if polylines are not equivalent.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def tojson(self) -> str:
        """
        Convert the polyline object to a JSON string.

        Returns:
            str: JSON representation of the polyline.

        Note:
            The 'aLine' and 'aPoint' attributes are excluded from JSON output
            to avoid circular references and reduce size.
        """
        aSkip = ["aLine", "aPoint"]

        obj = self.__dict__.copy()
        for sKey in aSkip:
            obj.pop(sKey, None)

        sJson = json.dumps(
            obj, sort_keys=True, indent=4, ensure_ascii=True, cls=PolylineClassEncoder
        )
        return sJson

    def towkt(self) -> str:
        """
        Convert the polyline object to a WKT string.

        Returns:
            str: WKT (Well-Known Text) representation of the polyline.
        """
        pGeometry = ogr.Geometry(ogr.wkbLineString)
        for i in range(self.nPoint):
            pGeometry.AddPoint(
                self.aPoint[i].dLongitude_degree, self.aPoint[i].dLatitude_degree
            )

        sWKT = pGeometry.ExportToWkt()
        pGeometry = None
        self.wkt = sWKT
        return sWKT

    def update_wkt(self) -> str:
        """
        Update and return the WKT string of the polyline.

        Returns:
            str: Updated WKT representation.
        """
        return self.towkt()
