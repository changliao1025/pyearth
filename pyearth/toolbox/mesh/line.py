import os
import json
from json import JSONEncoder
import numpy as np
import importlib.util
from typing import List, Tuple, Optional
from pyearth.toolbox.mesh.point import pypoint
from pyearth.toolbox.mesh.circle import pycircle
from pyearth.gis.geometry.calculate_intersect_on_great_circle import (
    calculate_intersect_on_great_circle,
)
from pyearth.gis.gdal.write.vector.gdal_export_point_to_vector_file import (
    export_point_as_polygon_file,
)
from pyearth.toolbox.mesh.algorithm.split_by_length import split_line_by_length
from pyearth.toolbox.mesh.algorithm.find_minimal_enclosing_polygon import (
    find_minimal_enclosing_polygon,
)
from pyearth.gis.geometry.calculate_angle_between_point import (
    calculate_angle_between_point,
)

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyearth.gis.geometry.kernel import calculate_distance_to_plane
else:
    from pyearth.gis.geometry.calculate_distance_to_plane import (
        calculate_distance_to_plane,
    )


# Constants for geometric calculations
ANGLE_THRESHOLD_COLLINEAR = (
    178.0  # Degrees - threshold for considering points collinear
)
ANGLE_THRESHOLD_PERPENDICULAR = 90.0  # Degrees - threshold for perpendicular check
DISTANCE_TOLERANCE = 1.0  # Meters - tolerance for point-on-line check
DISTANCE_PLANE_INITIAL = 9999.0  # Meters - initial large distance value


class LineClassEncoder(JSONEncoder):
    """Custom JSON encoder for pyline objects and their dependencies."""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float32):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pypoint):
            return json.loads(obj.tojson())
        return super().default(obj)


class pyline:
    """
    Represents a line segment (edge) between two geographic points.

    This class provides geometric operations for line segments including distance
    calculations, point-line relationships, line splitting, and buffer zone generation.
    All geographic calculations are performed on a spherical Earth model.

    Attributes:
        pPoint_start (pypoint): The starting point of the line.
        pPoint_end (pypoint): The ending point of the line.
        dLength (float): The geodesic length of the line in meters.
        pBound (Tuple[float, float, float, float]): Bounding box (lon_min, lat_min, lon_max, lat_max).
        lLineID (int): Optional identifier for the line (default: -1).
        lLineIndex (int): Optional index for the line (default: -1).

    Example:
        >>> from pyearth.toolbox.mesh.point import pypoint
        >>> p1 = pypoint({'dLongitude_degree': -122.4, 'dLatitude_degree': 37.8})
        >>> p2 = pypoint({'dLongitude_degree': -122.3, 'dLatitude_degree': 37.7})
        >>> line = pyline(p1, p2)
        >>> print(f"Line length: {line.dLength:.2f} meters")
    """

    def __init__(self, pPoint_start_in: pypoint, pPoint_end_in: pypoint):
        """
        Initialize a pyline object.

        Args:
            pPoint_start_in (pypoint): The starting point of the line.
            pPoint_end_in (pypoint): The ending point of the line.

        Raises:
            ValueError: If the two points are identical.
            TypeError: If inputs are not pypoint objects.
        """
        if not isinstance(pPoint_start_in, pypoint) or not isinstance(
            pPoint_end_in, pypoint
        ):
            raise TypeError("Both inputs must be pypoint objects.")

        if pPoint_start_in == pPoint_end_in:
            raise ValueError("Cannot create a line from two identical points.")

        self.pPoint_start: pypoint = pPoint_start_in
        self.pPoint_end: pypoint = pPoint_end_in
        self.dLength: float = self.calculate_length()
        self.pBound: Tuple[float, float, float, float] = self.calculate_line_bound()
        self.lLineID: int = -1
        self.lLineIndex: int = -1

    @classmethod
    def create(cls, pPoint_start_in: pypoint, pPoint_end_in: pypoint) -> "pyline":
        """
        Factory method to create a pyline object.

        Args:
            pPoint_start_in (pypoint): The starting point.
            pPoint_end_in (pypoint): The ending point.

        Returns:
            pyline: A pyline object.
        """
        return cls(pPoint_start_in, pPoint_end_in)

    def calculate_line_bound(self) -> Tuple[float, float, float, float]:
        """
        Calculate the bounding box of the line.

        Returns:
            Tuple[float, float, float, float]: (lon_min, lat_min, lon_max, lat_max) in degrees.
        """
        dLon_max = max(
            self.pPoint_start.dLongitude_degree, self.pPoint_end.dLongitude_degree
        )
        dLon_min = min(
            self.pPoint_start.dLongitude_degree, self.pPoint_end.dLongitude_degree
        )
        dLat_max = max(
            self.pPoint_start.dLatitude_degree, self.pPoint_end.dLatitude_degree
        )
        dLat_min = min(
            self.pPoint_start.dLatitude_degree, self.pPoint_end.dLatitude_degree
        )
        self.pBound = (dLon_min, dLat_min, dLon_max, dLat_max)
        return self.pBound

    def calculate_length(self) -> float:
        """
        Calculate the length of the line.

        Returns:
            float: The geodesic distance between start and end points in meters.
        """
        self.dLength = self.pPoint_start.calculate_distance(self.pPoint_end)
        return self.dLength

    def check_shared_point(self, other: "pyline") -> bool:
        """
        Check whether two lines share a common point.

        Args:
            other (pyline): Another line to compare with.

        Returns:
            bool: True if lines share at least one point, False otherwise.
        """
        return self.pPoint_start in (
            other.pPoint_start,
            other.pPoint_end,
        ) or self.pPoint_end in (other.pPoint_start, other.pPoint_end)

    def check_upstream(self, other: "pyline") -> bool:
        """
        Check if another line connects upstream to this line.

        Args:
            other (pyline): The potential upstream line.

        Returns:
            bool: True if other.pPoint_end equals this.pPoint_start.
        """
        return self.pPoint_start == other.pPoint_end

    def check_downstream(self, other: "pyline") -> bool:
        """
        Check if another line connects downstream from this line.

        Args:
            other (pyline): The potential downstream line.

        Returns:
            bool: True if this.pPoint_end equals other.pPoint_start.
        """
        return self.pPoint_end == other.pPoint_start

    def split_by_length(self, dLength_in: float) -> List["pyline"]:
        """
        Split the line into segments of specified maximum length.

        Args:
            dLength_in (float): Maximum length for each segment in meters.

        Returns:
            List[pyline]: List of line segments. Returns [self] if line is shorter than threshold.

        Raises:
            ValueError: If dLength_in is not positive.
        """
        if dLength_in <= 0:
            raise ValueError("Length threshold must be positive.")

        if self.dLength <= dLength_in:
            return [self]
        else:
            return split_line_by_length(self, dLength_in)

    def reverse(self) -> None:
        """
        Reverse the direction of the line by swapping start and end points.

        Note:
            This method recalculates dLength and pBound after reversing.
        """
        self.pPoint_start, self.pPoint_end = self.pPoint_end, self.pPoint_start
        # Recalculate since endpoints changed
        self.dLength = self.calculate_length()
        self.pBound = self.calculate_line_bound()

    def is_overlap(self, pLine_in: "pyline") -> bool:
        """
        Check if two lines represent the same segment (possibly reversed).

        Args:
            pLine_in (pyline): Another line to compare with.

        Returns:
            bool: True if lines have same endpoints (in any order), False otherwise.
        """
        return (
            self.pPoint_start == pLine_in.pPoint_start
            and self.pPoint_end == pLine_in.pPoint_end
        ) or (
            self.pPoint_start == pLine_in.pPoint_end
            and self.pPoint_end == pLine_in.pPoint_start
        )

    def check_point_on_line(self, pPoint_in: pypoint) -> Tuple[bool, float, float]:
        """
        Check if a point lies on this line segment.

        Args:
            pPoint_in (pypoint): The point to check.

        Returns:
            Tuple[bool, float, float]:
                - bool: True if point is on the line
                - float: Distance along line from start (-1 if not on line)
                - float: Perpendicular distance to line plane
        """
        dDistance = -1.0
        dDistance_plane = DISTANCE_PLANE_INITIAL
        on_edge = False

        if pPoint_in != self.pPoint_start and pPoint_in != self.pPoint_end:
            d1 = self.pPoint_start.calculate_distance(pPoint_in)
            d2 = self.pPoint_end.calculate_distance(pPoint_in)
            d3 = d1 + d2 - self.dLength
            angle3deg = calculate_angle_between_point(
                self.pPoint_start.dLongitude_degree,
                self.pPoint_start.dLatitude_degree,
                pPoint_in.dLongitude_degree,
                pPoint_in.dLatitude_degree,
                self.pPoint_end.dLongitude_degree,
                self.pPoint_end.dLatitude_degree,
            )

            dDistance_plane = calculate_distance_to_plane(
                self.pPoint_start.dLongitude_degree,
                self.pPoint_start.dLatitude_degree,
                pPoint_in.dLongitude_degree,
                pPoint_in.dLatitude_degree,
                self.pPoint_end.dLongitude_degree,
                self.pPoint_end.dLatitude_degree,
            )

            if angle3deg > ANGLE_THRESHOLD_COLLINEAR and d3 < DISTANCE_TOLERANCE:
                on_edge = True
                dDistance = d1

        return on_edge, dDistance, dDistance_plane

    def calculate_distance_to_point(self, pPoint_in: pypoint) -> Tuple[float, pypoint]:
        """
        Calculate the minimum distance from a point to this line segment.

        Args:
            pPoint_in (pypoint): The point to measure distance from.

        Returns:
            Tuple[float, pypoint]:
                - float: Minimum distance in meters
                - pypoint: The closest point on the line to pPoint_in
        """
        d1 = self.pPoint_start.calculate_distance(pPoint_in)
        d2 = self.pPoint_end.calculate_distance(pPoint_in)

        angle3deg = calculate_angle_between_point(
            self.pPoint_start.dLongitude_degree,
            self.pPoint_start.dLatitude_degree,
            pPoint_in.dLongitude_degree,
            pPoint_in.dLatitude_degree,
            self.pPoint_end.dLongitude_degree,
            self.pPoint_end.dLatitude_degree,
        )

        if angle3deg > ANGLE_THRESHOLD_PERPENDICULAR:
            dDistance_plane = calculate_distance_to_plane(
                self.pPoint_start.dLongitude_degree,
                self.pPoint_start.dLatitude_degree,
                pPoint_in.dLongitude_degree,
                pPoint_in.dLatitude_degree,
                self.pPoint_end.dLongitude_degree,
                self.pPoint_end.dLatitude_degree,
            )

            if dDistance_plane < d1 and dDistance_plane < d2:
                dLongitude_intersect, dLatitude_intersect = (
                    calculate_intersect_on_great_circle(
                        self.pPoint_start.dLongitude_degree,
                        self.pPoint_start.dLatitude_degree,
                        pPoint_in.dLongitude_degree,
                        pPoint_in.dLatitude_degree,
                        self.pPoint_end.dLongitude_degree,
                        self.pPoint_end.dLatitude_degree,
                    )
                )

                pPoint_out = pypoint(
                    {
                        "dLongitude_degree": dLongitude_intersect,
                        "dLatitude_degree": dLatitude_intersect,
                    }
                )
                dDistance_min = pPoint_out.calculate_distance(pPoint_in)

                if dDistance_min < d1 and dDistance_min < d2:
                    return dDistance_min, pPoint_out

        if d1 < d2:
            return d1, self.pPoint_start
        else:
            return d2, self.pPoint_end

    def calculate_angle_with(self, other: "pyline") -> float:
        """
        Calculate the angle between this line and another line.

        Uses the angle calculation function from the geometry module to compute
        the angle formed by the two line segments.

        Args:
            other (pyline): Another line to calculate angle with.

        Returns:
            float: Angle in degrees between the two lines.

        Raises:
            TypeError: If other is not a pyline object.

        Example:
            >>> line1 = pyline(p1, p2)
            >>> line2 = pyline(p2, p3)
            >>> angle = line1.calculate_angle_with(line2)
        """
        if not isinstance(other, pyline):
            raise TypeError(f"Expected pyline, got {type(other)}")

        return calculate_angle_between_point(
            self.pPoint_start.dLongitude_degree,
            self.pPoint_start.dLatitude_degree,
            self.pPoint_end.dLongitude_degree,
            self.pPoint_end.dLatitude_degree,
            other.pPoint_start.dLongitude_degree,
            other.pPoint_start.dLatitude_degree,
            other.pPoint_end.dLongitude_degree,
            other.pPoint_end.dLatitude_degree,
        )

    def calculate_buffer_zone_polygon(
        self,
        dRadius: float,
        nPoint: int = 36,
        sFilename_out: Optional[str] = None,
        sFolder_out: Optional[str] = None,
    ) -> Tuple[List[pypoint], List[pypoint], List[pypoint], List[pycircle]]:
        """
        Calculate a buffer zone polygon around this line.

        Args:
            dRadius (float): Buffer radius in meters.
            nPoint (int): Number of points to use for circular approximation (default: 36).
            sFilename_out (Optional[str]): Output file path for the buffer polygon.
            sFolder_out (Optional[str]): Output folder for intermediate buffer circles.

        Returns:
            Tuple containing:
                - List[pypoint]: Vertices of the buffer polygon
                - List[pypoint]: Center points used
                - List[pypoint]: All circle points generated
                - List[pycircle]: Circle objects created

        Raises:
            ValueError: If dRadius is not positive or nPoint is less than 3.
        """
        if dRadius <= 0:
            raise ValueError("Buffer radius must be positive.")
        if nPoint < 3:
            raise ValueError("Number of points must be at least 3.")

        if self.dLength < dRadius * 2.0:
            aEdge = [self]
        else:
            aEdge = self.split_by_length(dRadius * 2.0)

        aPoint_out = []
        aPoint_center = []
        aPoint_circle = []
        aCircle = []

        for i, pEdge in enumerate(aEdge):
            pPoint_start = pEdge.pPoint_start
            pPoint_end = pEdge.pPoint_end

            aPoint_center.extend([pPoint_start, pPoint_end])

            aPoint_start_buffer = pPoint_start.calculate_buffer_zone_circle(
                dRadius, nPoint
            )
            pEdge.pCircle_start = pycircle(pPoint_start, aPoint_start_buffer)
            aPoint_circle.extend(aPoint_start_buffer)

            aPoint_end_buffer = pPoint_end.calculate_buffer_zone_circle(dRadius, nPoint)
            pEdge.pCircle_end = pycircle(pPoint_end, aPoint_end_buffer)
            aPoint_circle.extend(aPoint_end_buffer)

            if sFolder_out:
                export_point_as_polygon_file(
                    aPoint_start_buffer,
                    os.path.join(sFolder_out, f"buffer_zone_start_{i}.geojson"),
                )
                export_point_as_polygon_file(
                    aPoint_end_buffer,
                    os.path.join(sFolder_out, f"buffer_zone_end_{i}.geojson"),
                )

            aCircle.extend([pEdge.pCircle_start, pEdge.pCircle_end])

        aLongitude_degree = [v.dLongitude_degree for v in aPoint_circle]
        aLatitude_degree = [v.dLatitude_degree for v in aPoint_circle]

        pPolygon_out = find_minimal_enclosing_polygon(
            aLongitude_degree, aLatitude_degree
        )

        aPoint_out = [
            pypoint({"dLongitude_degree": p[0], "dLatitude_degree": p[1]})
            for p in pPolygon_out
        ]

        if sFilename_out:
            export_point_as_polygon_file(aPoint_out, sFilename_out)

        return aPoint_out, aPoint_center, aPoint_circle, aCircle

    def __eq__(self, other: object) -> bool:
        """
        Check if two lines are equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if lines have identical start and end points in the same order.
        """
        if not isinstance(other, pyline):
            return NotImplemented
        return (
            self.pPoint_start == other.pPoint_start
            and self.pPoint_end == other.pPoint_end
        )

    def __ne__(self, other: object) -> bool:
        """
        Check if two lines are not equivalent.

        Args:
            other: Object to compare with.

        Returns:
            bool: True if lines are not equivalent.
        """
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self) -> int:
        """
        Generate hash for the line based on its endpoints.

        Returns:
            int: Hash value for the line.
        """
        return hash((self.pPoint_start, self.pPoint_end))

    def __repr__(self) -> str:
        """
        Return a detailed string representation for debugging.

        Returns:
            str: Developer-friendly representation.
        """
        return (
            f"pyline(start=({self.pPoint_start.dLongitude_degree:.6f}, "
            f"{self.pPoint_start.dLatitude_degree:.6f}), "
            f"end=({self.pPoint_end.dLongitude_degree:.6f}, "
            f"{self.pPoint_end.dLatitude_degree:.6f}), "
            f"length={self.dLength:.2f}m)"
        )

    def __str__(self) -> str:
        """
        Return a user-friendly string representation.

        Returns:
            str: Human-readable description.
        """
        return (
            f"Line from ({self.pPoint_start.dLongitude_degree:.4f}°, "
            f"{self.pPoint_start.dLatitude_degree:.4f}°) to "
            f"({self.pPoint_end.dLongitude_degree:.4f}°, "
            f"{self.pPoint_end.dLatitude_degree:.4f}°)"
        )

    def tojson(self) -> str:
        """
        Convert the line object to a JSON string.

        Returns:
            str: JSON representation of the line.
        """
        return json.dumps(self, cls=LineClassEncoder, sort_keys=True, indent=4)
