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

def split_flowline_by_length(aFlowline_in, dDistance):

    aFlowline_out=list()
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        pFlowline_out = pFlowline.split_by_length(dDistance)
        aFlowline_out.append(pFlowline_out)

    return aFlowline_out

def split_edge_by_length(pEdge_in, dLength_in, tolerance=1e-6):
    """
    Split an edge into smaller segments with maximum length constraint.

    Uses optimal segmentation instead of recursive binary division for better
    performance and more even segment distribution.

    Args:
        pEdge_in: Input edge to split
        dLength_in: Maximum length for each segment
        tolerance: Relative tolerance for length comparison

    Returns:
        List of edge segments

    Raises:
        ValueError: If length threshold is not positive
        AttributeError: If edge doesn't have required attributes
    """
    from pyflowline.classes.edge import pyedge

    # Input validation
    if dLength_in <= 0:
        raise ValueError("Length threshold must be positive")

    if not hasattr(pEdge_in, 'dLength'):
        raise AttributeError("Edge must have dLength attribute")

    dLength_total = pEdge_in.dLength

    # Early return for edges that are already short enough
    if dLength_total <= dLength_in * (1 + tolerance):
        return [pEdge_in]

    # Calculate optimal number of segments
    nSegments = int(np.ceil(dLength_total / dLength_in))

    pVertex_start = pEdge_in.pVertex_start
    pVertex_end = pEdge_in.pVertex_end

    n1 = pVertex_start.toNvector()
    n2 = pVertex_end.toNvector()

    aEdge_out = []

    # Create evenly spaced segments using spherical interpolation
    for i in range(nSegments):
        t1 = i / nSegments
        t2 = (i + 1) / nSegments

        # Handle endpoints to avoid unnecessary interpolation
        if i == 0:
            pVertex_start_seg = pVertex_start
        else:
            mid1 = slerp(n1, n2, t1)
            pVertex_start_seg = mid1.toLatLon()

        if i == nSegments - 1:
            pVertex_end_seg = pVertex_end
        else:
            mid2 = slerp(n1, n2, t2)
            pVertex_end_seg = mid2.toLatLon()

        # Create edge segment
        pEdge = pyedge(pVertex_start_seg, pVertex_end_seg)
        aEdge_out.append(pEdge)

    return aEdge_out

