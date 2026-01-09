import numpy as np


def cross_product_2d(p_prev, p_curr, p_next):
    """
    Calculates the 2D cross product of vectors defined by three consecutive points.

    The sign of the result indicates the direction of the turn at p_curr:
    > 0: Left turn (Counter-Clockwise, CCW)
    < 0: Right turn (Clockwise, CW)
    = 0: Collinear (straight line)

    Vectors: A = P_curr - P_prev, B = P_next - P_curr
    """
    # Vector A components (Incoming edge)
    ax = p_curr[0] - p_prev[0]
    ay = p_curr[1] - p_prev[1]

    # Vector B components (Outgoing edge)
    bx = p_next[0] - p_curr[0]
    by = p_next[1] - p_curr[1]

    # 2D Cross Product (Ax * By) - (Ay * Bx)
    return (ax * by) - (ay * bx)


def is_convex_polygon(polygon_vertices):
    """
    Checks if the polygon defined by a list of vertices is convex.

    Args:
        polygon_vertices (list of tuples or list of lists): List of [x, y] coordinates.

    Returns:
        bool: True if the polygon is convex, False otherwise.
    """
    n = len(polygon_vertices)

    # Polygons must have at least 3 vertices to be considered (or 4 if the start/end point is repeated)
    if n < 3:
        return True

    # We will store the sign of the first non-zero cross product found (i.e., the expected turn direction)
    # 0: Undetermined (only collinear points found so far)
    # 1: Expected CCW (Positive)
    # -1: Expected CW (Negative)
    first_sign = 0

    # Iterate through all vertices (i) as the middle point,
    # checking the turn formed by (i-1), i, and (i+1).
    for i in range(n):

        # Use modulo (%) for wrapping around the list (P_n is P_0, P_0 is P_n)
        p_prev = polygon_vertices[(i - 1 + n) % n]
        p_curr = polygon_vertices[i]
        p_next = polygon_vertices[(i + 1) % n]

        # Calculate the direction of the turn at p_curr
        cp = cross_product_2d(p_prev, p_curr, p_next)

        # 4. Check for consistency
        if cp != 0:
            current_sign = 1 if cp > 0 else -1

            if first_sign == 0:
                # Set the required sign based on the first definite turn
                first_sign = current_sign
            elif current_sign != first_sign:
                # If the sign flips, the polygon is concave because it has an internal angle > 180 degrees.
                return False

    # If the loop finishes without the sign flipping, the polygon is convex.
    return True
