import numpy as np


def slerp(n1, n2, t):
    """
    Spherical linear interpolation between two n-vectors.

    Args:
        n1: First n-vector
        n2: Second n-vector
        t: Interpolation parameter (0 to 1)

    Returns:
        Interpolated n-vector
    """
    dot_product = n1.dot(n2)

    # Handle nearly parallel vectors
    if abs(dot_product) > 0.9995:
        # Linear interpolation for very close points
        result = n1 * (1 - t) + n2 * t
        return result.normalize()

    # Standard slerp
    theta = np.arccos(np.clip(dot_product, -1, 1))
    sin_theta = np.sin(theta)

    factor1 = np.sin((1 - t) * theta) / sin_theta
    factor2 = np.sin(t * theta) / sin_theta

    return n1 * factor1 + n2 * factor2


def split_polyline_by_length(aFlowline_in, dDistance):

    aPolyline_out = list()
    nPolyline = len(aFlowline_in)
    for i in range(nPolyline):
        pPolyline = aFlowline_in[i]
        pPolyline_out = pPolyline.split_by_length(dDistance)
        aPolyline_out.append(pPolyline_out)

    return aPolyline_out


def split_line_by_length(pLine_in, dLength_in, tolerance=1e-6):
    """
    Split a line into smaller segments with maximum length constraint.

    Uses optimal segmentation instead of recursive binary division for better
    performance and more even segment distribution.

    Args:
        pLine_in: Input line to split
        dLength_in: Maximum length for each segment
        tolerance: Relative tolerance for length comparison

    Returns:
        List of line segments

    Raises:
        ValueError: If length threshold is not positive
        AttributeError: If line doesn't have required attributes
    """
    from pyearth.toolbox.mesh.line import pyline

    # Input validation
    if dLength_in <= 0:
        raise ValueError("Length threshold must be positive")

    if not hasattr(pLine_in, "dLength"):
        raise AttributeError("Line must have dLength attribute")

    dLength_total = pLine_in.dLength

    # Early return for lines that are already short enough
    if dLength_total <= dLength_in * (1 + tolerance):
        return [pLine_in]

    # Calculate optimal number of segments
    nSegments = int(np.ceil(dLength_total / dLength_in))

    pPoint_start = pLine_in.pPoint_start
    pPoint_end = pLine_in.pPoint_end

    n1 = pPoint_start.toNvector()
    n2 = pPoint_end.toNvector()

    aLine_out = []

    # Create evenly spaced segments using spherical interpolation

    for i in range(nSegments):
        t1 = i / nSegments
        t2 = (i + 1) / nSegments

        # Handle endpoints to avoid unnecessary interpolation
        if i == 0:
            pPoint_start_seg = pPoint_start
        else:
            mid1 = slerp(n1, n2, t1)
            pPoint_start_seg = mid1.toLatLon()

        if i == nSegments - 1:
            pPoint_end_seg = pPoint_end
        else:
            mid2 = slerp(n1, n2, t2)
            pPoint_end_seg = mid2.toLatLon()

        # Check for degenerate segments (same start and end points)
        if pPoint_start_seg == pPoint_end_seg:
            # Skip degenerate segments or use a small offset
            print(f"Warning: Segment {i} has identical start and end points")
            print(
                f"  Start: ({pPoint_start_seg.dLongitude_degree}, {pPoint_start_seg.dLatitude_degree})"
            )
            print(
                f"  End: ({pPoint_end_seg.dLongitude_degree}, {pPoint_end_seg.dLatitude_degree})"
            )
            print(f"  t1={t1}, t2={t2}, nSegments={nSegments}")
            continue

        # Create line segment
        pLine = pyline(pPoint_start_seg, pPoint_end_seg)
        aLine_out.append(pLine)

    return aLine_out
