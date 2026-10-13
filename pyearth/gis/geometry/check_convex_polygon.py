"""Check if a 2D polygon is convex."""

from typing import List, Sequence, Tuple, Union
import numpy as np

__all__ = [
    "is_convex_polygon",
    "cross_product_2d",
]


def cross_product_2d(p_prev: Sequence[float], p_curr: Sequence[float], p_next: Sequence[float]) -> float:
    """Calculate the 2D cross product (turn direction) of vectors defined by three consecutive points.

    Returns:
    > 0: Left turn (Counter-Clockwise, CCW)
    < 0: Right turn (Clockwise, CW)
    = 0: Collinear
    """
    ax = p_curr[0] - p_prev[0]
    ay = p_curr[1] - p_prev[1]
    bx = p_next[0] - p_curr[0]
    by = p_next[1] - p_curr[1]
    return float((ax * by) - (ay * bx))


def is_convex_polygon(polygon_vertices: Sequence[Sequence[float]]) -> bool:
    """Check if the polygon defined by a sequence of vertices is convex.

    Parameters
    ----------
    polygon_vertices : Sequence[Sequence[float]]
        List or array of [x, y] coordinates. Can be open or closed.

    Returns
    -------
    bool
        True if the polygon is convex, False otherwise.
    """
    vertices = list(polygon_vertices)
    if len(vertices) > 1 and np.allclose(vertices[0], vertices[-1]):
        vertices = vertices[:-1]

    n = len(vertices)
    if n < 3:
        return True

    first_sign = 0

    for i in range(n):
        p_prev = vertices[(i - 1 + n) % n]
        p_curr = vertices[i]
        p_next = vertices[(i + 1) % n]

        cp = cross_product_2d(p_prev, p_curr, p_next)

        if not np.isclose(cp, 0.0, atol=1e-12):
            current_sign = 1 if cp > 0 else -1

            if first_sign == 0:
                first_sign = current_sign
            elif current_sign != first_sign:
                return False

    return True

