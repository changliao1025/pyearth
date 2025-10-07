import math
from typing import List, Tuple

# Define the coordinate type alias for clarity (Lon, Lat)
Coord = Tuple[float, float]

def reorder_idl_polygon(vertices: List[Coord]) -> List[Coord]:
    """
    Reorders the vertices of an International Date Line (IDL)-crossing polygon
    to prevent self-intersection.

    This strategy assumes the polygon is tightly clustered around the 180/-180 meridian,
    and sorts the vertices to trace the perimeter consistently: up the Eastern side,
    across, down the Western side, and back across.

    Args:
        vertices: A list of (Longitude, Latitude) tuples. The list is expected
                  to be closed (first point equals the last point).

    Returns:
        A new list of (Longitude, Latitude) tuples forming a valid, closed polygon.
    """
    if not vertices:
        return []

    # 1. Strip the closing vertex for processing
    # We will re-add it correctly at the end.
    unique_vertices = vertices[:-1]

    # 2. Separate points into East (positive Lon) and West (negative Lon)
    east_points: List[Coord] = []
    west_points: List[Coord] = []

    # Note: Longitude 0.0 (Prime Meridian) is included in the East group.
    for lon, lat in unique_vertices:
        if lon >= 0:
            east_points.append((lon, lat))
        else:
            west_points.append((lon, lat))

    # Handle edge case where all points are on one side (no IDL crossing)
    if not east_points or not west_points:
        # If it doesn't cross, it's not at risk of IDL self-intersection.
        return vertices

    # 3. Sort the points to define a perimeter tracing path

    # Sort East points by Latitude (ascending) -> Traces East boundary from South to North
    east_points.sort(key=lambda p: p[1])

    # Sort West points by Latitude (descending) -> Traces West boundary from North to South
    west_points.sort(key=lambda p: p[1], reverse=True)

    # 4. Construct the new, valid linear ring
    # Path: [South to North on East] -> [North to South on West]
    new_vertices = east_points + west_points

    # 5. Close the polygon by appending the starting point
    new_vertices.append(new_vertices[0])

    #convert it back to ccw
    new_vertices = new_vertices[::-1]

    return new_vertices