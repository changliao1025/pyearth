import numpy as np
import logging

# Configure logger for this module
logger = logging.getLogger(__name__)


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
    from pyearth.toolbox.mesh.nvector import pynvector

    use_high_precision = getattr(n1, "use_high_precision", False)
    dtype = np.float128 if use_high_precision else np.float64
    eps = dtype(1e-18 if use_high_precision else 1e-12)

    dot_product = dtype(np.clip(n1.dot(n2), -1.0, 1.0))

    # Near-parallel vectors: lerp is stable and avoids division by tiny sin(theta).
    if dot_product > dtype(1.0) - eps:
        result = n1 * (dtype(1.0) - dtype(t)) + n2 * dtype(t)
        return result.normalize()

    # Near-antiparallel vectors: avoid lerp because midpoint can collapse to zero.
    if dot_product < -dtype(1.0) + eps:
        abs_x = abs(n1.dX)
        abs_y = abs(n1.dY)
        abs_z = abs(n1.dZ)

        # Build a deterministic orthogonal direction to n1.
        if abs_x <= abs_y and abs_x <= abs_z:
            ortho = pynvector(
                {"x": dtype(0.0), "y": -n1.dZ, "z": n1.dY},
                use_high_precision=use_high_precision,
            )
        elif abs_y <= abs_x and abs_y <= abs_z:
            ortho = pynvector(
                {"x": -n1.dZ, "y": dtype(0.0), "z": n1.dX},
                use_high_precision=use_high_precision,
            )
        else:
            ortho = pynvector(
                {"x": -n1.dY, "y": n1.dX, "z": dtype(0.0)},
                use_high_precision=use_high_precision,
            )

        ortho = ortho.normalize()
        theta = dtype(np.pi) * dtype(t)
        return (n1 * np.cos(theta) + ortho * np.sin(theta)).normalize()

    # Standard slerp
    theta = np.arccos(dot_product)
    sin_theta = np.sin(theta)

    if abs(sin_theta) <= eps:
        result = n1 * (dtype(1.0) - dtype(t)) + n2 * dtype(t)
        return result.normalize()

    factor1 = np.sin((1 - t) * theta) / sin_theta
    factor2 = np.sin(t * theta) / sin_theta

    return (n1 * factor1 + n2 * factor2).normalize()


def split_polyline_by_length(aFlowline_in, dDistance):

    aPolyline_out = list()
    nPolyline = len(aFlowline_in)
    for i in range(nPolyline):
        pPolyline = aFlowline_in[i]
        pPolyline_out = pPolyline.split_by_length(dDistance)
        aPolyline_out.append(pPolyline_out)

    return aPolyline_out


def split_line_by_length(pLine_in, dLength_in, tolerance=1e-6, use_high_precision=True):
    """
    Split a line into smaller segments with maximum length constraint.

    Uses optimal segmentation with spherical interpolation. Supports high-precision
    mode (float128) to reduce numerical errors in the coordinate conversion chain.

    Args:
        pLine_in: Input line to split
        dLength_in: Maximum length for each segment (meters)
        tolerance: Relative tolerance for length comparison (default: 1e-6)
        use_high_precision: Use float128 for conversions (default: True)

    Returns:
        List of line segments

    Raises:
        ValueError: If length threshold is not positive
        AttributeError: If line doesn't have required attributes

    Note:
        When use_high_precision=True, the conversion chain uses numpy.float128
        to maintain precision during spherical interpolation. This significantly
        reduces degenerate segment warnings but has a small performance cost (~10-30%).

        The adaptive segmentation limits the number of segments based on numerical
        precision to prevent creating segments that cannot be distinguished.
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

    # Adaptive segmentation: limit based on numerical precision
    # Estimate the angular span on the sphere
    dRadius_earth = 6371000.0  # meters (mean Earth radius)
    dAngle_total = dLength_total / dRadius_earth  # approximate angle in radians

    # Minimum distinguishable angle based on precision
    # float128: ~33 decimal digits, can distinguish ~1e-30 radians
    # float64: ~15 decimal digits, can distinguish ~1e-15 radians
    # Use conservative estimates to ensure points remain distinguishable
    dAngle_min = 1e-14 if use_high_precision else 1e-12
    nSegments_max = max(1, int(dAngle_total / dAngle_min))

    if nSegments > nSegments_max:
        logger.info(
            f"Reducing segments from {nSegments} to {nSegments_max} "
            f"due to precision limits (line length={dLength_total:.2f}m, "
            f"angle={np.degrees(dAngle_total):.6e}°)"
        )
        nSegments = nSegments_max

    pPoint_start = pLine_in.pPoint_start
    pPoint_end = pLine_in.pPoint_end

    # Use high precision for n-vector conversion if requested
    n1 = pPoint_start.toNvector(use_high_precision=use_high_precision)
    n2 = pPoint_end.toNvector(use_high_precision=use_high_precision)

    aLine_out = []

    # Build segments from the last accepted endpoint so skipped degenerate
    # interpolation intervals are merged into the next valid segment.
    dThreshold_segment_m = 1e-3  # 1 mm
    pPoint_start_seg = pPoint_start

    for i in range(nSegments):
        t1 = i / nSegments
        t2 = (i + 1) / nSegments

        if i == nSegments - 1:
            pPoint_end_seg = pPoint_end
        else:
            mid2 = slerp(n1, n2, t2)
            pPoint_end_seg = mid2.toLatLon()

        # Check for degenerate segments (nearly same start and end points)
        dSeg_candidate = pPoint_start_seg.calculate_distance(pPoint_end_seg)
        if dSeg_candidate <= dThreshold_segment_m:
            logger.debug(
                f"Skipping degenerate interval {i}/{nSegments}: "
                f"coords=({pPoint_start_seg.dLongitude_degree:.15f}, "
                f"{pPoint_start_seg.dLatitude_degree:.15f}), "
                f"distance={dSeg_candidate:.6e}m, "
                f"t1={t1:.6f}, t2={t2:.6f}, "
                f"precision={'float128' if use_high_precision else 'float64'}"
            )
            print(
                f"Skipping degenerate interval {i}/{nSegments}: "
                f"coords=({pPoint_start_seg.dLongitude_degree:.15f}, "
                f"{pPoint_start_seg.dLatitude_degree:.15f}), "
                f"distance={dSeg_candidate:.6e}m, "
                f"t1={t1:.6f}, t2={t2:.6f}, "
                f"precision={'float128' if use_high_precision else 'float64'}"
            )
            continue

        # Create line segment with error handling
        try:
            pLine = pyline(pPoint_start_seg, pPoint_end_seg)
            aLine_out.append(pLine)
            pPoint_start_seg = pPoint_end_seg
        except ValueError as e:
            logger.debug(f"Skipping interval {i}/{nSegments}: {str(e)}")
            print(f"Skipping interval {i}/{nSegments}: {str(e)}")
            continue

    return aLine_out
