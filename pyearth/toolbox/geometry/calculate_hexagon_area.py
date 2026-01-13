"""
Regular Hexagon Area Calculation

This module provides functionality for calculating the area of a regular hexagon
given its edge length. Regular hexagons are six-sided polygons with all sides of
equal length and all interior angles equal to 120 degrees.

Main Functions
--------------
calculate_hexagon_area : Calculate area of regular hexagon from edge length

Key Features
------------
- Accurate area calculation using standard geometric formula
- Supports scalar and array inputs via NumPy
- Input validation with informative error messages
- Type hints for better code clarity
- Professional logging for debugging
- Handles edge cases (zero, negative values)

Use Cases
---------
1. **Hexagonal Grids**: Calculate cell areas in hexagonal grid systems
2. **Geographic Modeling**: Compute areas for H3 or DGGS hexagonal cells
3. **Game Development**: Calculate hexagonal tile areas
4. **Spatial Analysis**: Determine coverage areas in hexagonal tessellations
5. **Materials Science**: Calculate cross-sectional areas of hexagonal structures
6. **Urban Planning**: Analyze hexagonal city block patterns
7. **Data Visualization**: Size hexagonal bins in hex bin plots

Technical Details
-----------------
For a regular hexagon with edge length s, the area is calculated using:

    Area = (3√3/2) × s²

This can also be expressed as:
    Area ≈ 2.598076211353316 × s²

Derivation:
- A regular hexagon can be divided into 6 equilateral triangles
- Each triangle has side length s and area (√3/4) × s²
- Total area = 6 × (√3/4) × s² = (6√3/4) × s² = (3√3/2) × s²

Alternative formulas (equivalent):
- From apothem (a): Area = 2√3 × a²
- From circumradius (R): Area = (3√3/2) × R²
- From inradius/apothem (a): Area = 2√3 × a²

Geometric Properties:
- Interior angle: 120°
- Exterior angle: 60°
- Number of sides: 6
- Number of vertices: 6
- Perimeter: 6s
- Apothem (inradius): a = (√3/2) × s
- Circumradius: R = s

Performance Characteristics
---------------------------
- Time Complexity: O(1) for scalar, O(N) for arrays
- Space Complexity: O(1) for scalar, O(N) for arrays
- Vectorized operations for array inputs

Dependencies
------------
- numpy: Numerical operations and array support
- logging: For error reporting and debugging

See Also
--------
- numpy.sqrt: Square root function
- calculate_polygon_area: General polygon area calculation
"""

import numpy as np
import logging
from typing import Union

# Configure logging
logger = logging.getLogger(__name__)


def calculate_hexagon_area(
    dLength_edge_in: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """
    Calculate the area of a regular hexagon from its edge length.

    This function computes the area of a regular hexagon (all sides equal, all
    angles 120°) using the standard geometric formula. The calculation works
    for both scalar values and NumPy arrays.

    Parameters
    ----------
    dLength_edge_in : float or np.ndarray
        The length of the hexagon edge (side). Must be non-negative.
        Can be a scalar float or a NumPy array of floats.

        - For scalar: Returns scalar area
        - For array: Returns array of areas (vectorized operation)

    Returns
    -------
    float or np.ndarray
        The calculated hexagon area. Return type matches input type.

        - For edge length s: Area = (3√3/2) × s² ≈ 2.598076 × s²
        - For zero edge length: Returns 0.0
        - For array input: Returns array of same shape

    Raises
    ------
    TypeError
        If dLength_edge_in is None or not a numeric type.
    ValueError
        If dLength_edge_in contains negative values.

    Notes
    -----
    1. **Formula**: Area = (3√3/2) × s² where s is the edge length
    2. **Geometric Basis**: A regular hexagon can be divided into 6 equilateral
       triangles, each with area (√3/4) × s²
    3. **Numerical Constant**: 3√3/2 ≈ 2.598076211353316
    4. **Zero Handling**: Edge length of 0 correctly returns area of 0
    5. **Array Support**: Fully vectorized for NumPy arrays (no loops)
    6. **Unit Consistency**: Output units are square of input units (e.g., m² for m)
    7. **Precision**: Uses NumPy's floating-point precision (typically float64)
    8. **Alternative Metrics**:
       - From apothem a: Area = 2√3 × a² (where a = (√3/2) × s)
       - From circumradius R: Area = (3√3/2) × R² (where R = s)
    9. **Related Formulas**:
       - Perimeter: P = 6s
       - Apothem (inradius): a = s × √3/2 ≈ 0.866s
       - Circumradius: R = s
    10. **Common Applications**: Hexagonal grids (H3, DGGS), game boards,
        spatial tessellations, material science cross-sections

    Examples
    --------
    Calculate area of hexagon with edge length 1 meter:

    >>> area = calculate_hexagon_area(1.0)
    >>> print(f"Area: {area:.6f} m²")
    Area: 2.598076 m²

    Calculate area of hexagon with edge length 10 kilometers:

    >>> area_km = calculate_hexagon_area(10.0)
    >>> print(f"Area: {area_km:.2f} km²")
    Area: 259.81 km²

    Calculate areas for multiple hexagons (H3 cell analysis):

    >>> edge_lengths = np.array([1.0, 2.0, 5.0, 10.0])
    >>> areas = calculate_hexagon_area(edge_lengths)
    >>> for edge, area in zip(edge_lengths, areas):
    ...     print(f"Edge {edge:4.1f} -> Area {area:8.3f}")
    Edge  1.0 -> Area    2.598
    Edge  2.0 -> Area   10.392
    Edge  5.0 -> Area   64.952
    Edge 10.0 -> Area  259.808

    Verify zero edge length:

    >>> area_zero = calculate_hexagon_area(0.0)
    >>> print(f"Area of zero edge: {area_zero}")
    Area of zero edge: 0.0

    Calculate area with different units (inches to square inches):

    >>> edge_inches = 6.0
    >>> area_sq_inches = calculate_hexagon_area(edge_inches)
    >>> print(f"{edge_inches} inch edge → {area_sq_inches:.2f} in²")
    6.0 inch edge → 93.53 in²

    Array calculation for hexagonal grid:

    >>> # Create 3x3 grid of edge lengths
    >>> grid_edges = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=float)
    >>> grid_areas = calculate_hexagon_area(grid_edges)
    >>> print("Grid areas:")
    >>> print(grid_areas)
    Grid areas:
    [[  2.598  10.392  23.383]
     [ 41.569  64.952  93.531]
     [127.305 166.277 210.444]]

    See Also
    --------
    numpy.sqrt : Square root function used in calculation

    References
    ----------
    .. [1] Weisstein, Eric W. "Regular Hexagon." From MathWorld--A Wolfram Web
           Resource. https://mathworld.wolfram.com/RegularHexagon.html
    .. [2] H3: Uber's Hexagonal Hierarchical Spatial Index
           https://h3geo.org/
    """
    # Validate input is not None
    if dLength_edge_in is None:
        error_msg = "Edge length cannot be None"
        logger.error(error_msg)
        raise TypeError(error_msg)

    # Convert to numpy array for consistent handling
    edge_array = np.asarray(dLength_edge_in)

    # Validate numeric type
    if not np.issubdtype(edge_array.dtype, np.number):
        error_msg = f"Edge length must be numeric, got {type(dLength_edge_in).__name__}"
        logger.error(error_msg)
        raise TypeError(error_msg)

    # Validate non-negative values
    if np.any(edge_array < 0):
        error_msg = "Edge length must be non-negative"
        logger.error(error_msg)
        raise ValueError(error_msg)

    # Calculate area using correct formula: Area = (3√3/2) × s²
    # This is equivalent to: Area = 1.5 × √3 × s²
    area = 1.5 * np.sqrt(3.0) * edge_array**2

    # Return scalar if input was scalar
    if np.isscalar(dLength_edge_in):
        return float(area)

    return area
