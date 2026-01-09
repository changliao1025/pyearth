"""
Convert longitude values between different coordinate range conventions.

This module provides utilities for converting longitude values between two common
coordinate systems used in geographic data:
    - [-180°, 180°] range: Standard geographic coordinates (western/eastern hemispheres)
    - [0°, 360°] range: Alternative representation (continuous eastward from prime meridian)

The conversions use modulo arithmetic to handle values outside the target range
and ensure consistent wrapping behavior.

Functions
---------
convert_360_to_180
    Convert longitude value(s) from [0°, 360°] to [-180°, 180°] range.
    Handles both scalar values and arrays automatically.
convert_180_to_360
    Convert a single longitude value from [-180°, 180°] to [0°, 360°] range.

Notes
-----
**Coordinate System Conventions**:

1. **[-180°, 180°] System** (Standard):
   - Western hemisphere: negative values (-180° to 0°)
   - Eastern hemisphere: positive values (0° to 180°)
   - Prime meridian: 0°
   - International Date Line: ±180° (discontinuity)
   - Used by: Most GIS software, GPS devices, standard projections

2. **[0°, 360°] System** (Alternative):
   - Continuous eastward from prime meridian
   - Prime meridian: 0° (or 360°)
   - International Date Line: 180°
   - No negative values, no discontinuity at ±180°
   - Used by: Some climate models, oceanographic data, global grids

**Conversion Formula**:

The conversion uses modulo arithmetic to ensure consistent wrapping:
    - 360° to 180°: ((lon + 180) % 360) - 180
    - 180° to 360°: (lon + 360) % 360

**Edge Cases**:

- Values at boundaries (0°, ±180°, 360°) are handled consistently
- Values outside the expected range are automatically wrapped
- The 180° meridian can be represented as either 180° or -180° in the
  [-180°, 180°] system; this implementation uses -180°

Examples
--------
Convert individual longitude values:

>>> lon_360 = 270.0  # 270°E in [0, 360] system
>>> lon_180 = convert_360_to_180(lon_360)
>>> print(lon_180)
-90.0  # 90°W in [-180, 180] system

>>> lon_180 = -120.0  # 120°W
>>> lon_360 = convert_180_to_360(lon_180)
>>> print(lon_360)
240.0  # 240°E

Convert arrays of longitudes:

>>> import numpy as np
>>> lons_360 = np.array([0, 90, 180, 270, 360])
>>> lons_180 = convert_360_to_180(lons_360)
>>> print(lons_180)
[0, 90, -180, -90, 0]

See Also
--------
numpy.mod : Modulo operation for arrays
pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon : Handle IDL-crossing polygons

References
----------
.. [1] IDL Coyote longitude conversion tips:
       http://www.idlcoyote.com/map_tips/lonconvert.html
"""

from typing import Union, overload
import numpy as np


@overload
def convert_360_to_180(longitude_degree_in: float) -> float: ...


@overload
def convert_360_to_180(
    longitude_degree_in: Union[np.ndarray, list, tuple],
) -> np.ndarray: ...


def convert_360_to_180(
    longitude_degree_in: Union[float, np.ndarray, list, tuple],
) -> Union[float, np.ndarray]:
    """
    Convert longitude from [0°, 360°] range to [-180°, 180°] range.

    This function converts longitude values from the continuous eastward system
    [0°, 360°] to the standard western/eastern hemisphere system [-180°, 180°].
    The conversion uses modulo arithmetic to handle values outside the expected
    range and ensure consistent wrapping.

    The function automatically handles both scalar values and arrays, returning
    the same type as the input.

    Parameters
    ----------
    longitude_degree_in : float or array-like
        Input longitude in degrees. Can be:
        - A single float value
        - A numpy array
        - A list or tuple of values

        Values are expected to be in the range [0°, 360°], but values outside
        this range are automatically wrapped. Can be any real numbers.

    Returns
    -------
    float or numpy.ndarray
        Longitude in degrees in the range [-180°, 180°].
        - If input is a scalar float, returns a float
        - If input is an array, list, or tuple, returns a numpy array

        The value -180° and 180° represent the same meridian (International Date Line).

    Notes
    -----
    **Conversion Formula**:

        longitude_out = ((longitude_in + 180) % 360) - 180

    This formula works by:
    1. Shifting input by +180° to [180°, 540°] range
    2. Taking modulo 360 to wrap to [0°, 360°] range
    3. Shifting back by -180° to [-180°, 180°] range

    **Examples of Conversion**:

    - 0° → 0° (Prime Meridian)
    - 90° → 90° (Eastern hemisphere)
    - 180° → -180° (International Date Line)
    - 270° → -90° (Western hemisphere)
    - 360° → 0° (wraps around to Prime Meridian)

    **Handling Out-of-Range Values**:

    Values outside [0°, 360°] are automatically wrapped:
    - 370° → 10° (wraps around once)
    - -10° → -10° (already in target range)
    - 720° → 0° (wraps around twice)

    **Array Operations**:

    For array inputs, the function uses vectorized numpy operations which are
    typically 10-100x faster than Python loops due to C-level optimization.
    Arrays, lists, and tuples are all handled efficiently.

    Warnings
    --------
    - For longitude values at exactly 180°, this function returns -180° rather
      than 180°. Both values represent the same meridian (International Date Line),
      but the choice affects how geometries crossing the IDL are represented.
    - Array inputs are converted to numpy arrays, which may create a copy
    - NaN and Inf values in arrays are preserved in the output

    Examples
    --------
    Convert single scalar values:

    >>> convert_360_to_180(0.0)
    0.0
    >>> convert_360_to_180(90.0)
    90.0
    >>> convert_360_to_180(180.0)
    -180.0
    >>> convert_360_to_180(270.0)
    -90.0
    >>> convert_360_to_180(360.0)
    0.0

    Handle values outside expected range:

    >>> convert_360_to_180(450.0)
    90.0
    >>> convert_360_to_180(-30.0)
    -30.0

    Convert numpy arrays:

    >>> import numpy as np
    >>> lons = np.array([0, 90, 180, 270, 360])
    >>> result = convert_360_to_180(lons)
    >>> print(result)
    [0, 90, -180, -90, 0]

    Convert from a list:

    >>> lons_list = [45.0, 135.0, 225.0, 315.0]
    >>> result = convert_360_to_180(lons_list)
    >>> print(result)
    [45.0, 135.0, -135.0, -45.0]

    Convert global grid longitudes:

    >>> # Create a global grid from 0° to 360° with 5° spacing
    >>> lons_360 = np.arange(0, 361, 5)
    >>> lons_180 = convert_360_to_180(lons_360)
    >>> print(lons_180[:5])  # First 5 values
    [0, 5, 10, 15, 20]
    >>> print(lons_180[-5:])  # Last 5 values
    [-20, -15, -10, -5, 0]

    Preserve special values in arrays:

    >>> lons = np.array([90.0, np.nan, 270.0, np.inf])
    >>> result = convert_360_to_180(lons)
    >>> print(result)
    [90.0, nan, -90.0, inf]

    Convert 2D arrays (e.g., meshgrid):

    >>> lon_2d = np.array([[0, 180], [90, 270]])
    >>> result = convert_360_to_180(lon_2d)
    >>> print(result)
    [[0, -180], [90, -90]]

    Convert coordinates from global dataset:

    >>> # Ocean point in Pacific at 200°E
    >>> lon = convert_360_to_180(200.0)
    >>> print(f"{lon}°")  # -160°W
    -160.0°

    See Also
    --------
    convert_180_to_360 : Inverse conversion ([-180, 180] to [0, 360])
    numpy.asarray : Convert input to array
    numpy.mod : Modulo operation for arrays

    References
    ----------
    .. [1] IDL Coyote longitude conversion:
           http://www.idlcoyote.com/map_tips/lonconvert.html
    """
    # Check if input is a scalar (float or int)
    if isinstance(longitude_degree_in, (int, float)):
        longitude_out = ((longitude_degree_in + 180) % 360) - 180
        return float(longitude_out)

    # Handle array-like inputs (numpy arrays, lists, tuples)
    longitude_array = np.asarray(longitude_degree_in)
    longitude_out = ((longitude_array + 180) % 360) - 180
    return longitude_out


def convert_180_to_360(longitude_degree_in: float) -> float:
    """
    Convert longitude from [-180°, 180°] range to [0°, 360°] range.

    This function converts longitude values from the standard western/eastern
    hemisphere system [-180°, 180°] to the continuous eastward system [0°, 360°].
    The conversion uses modulo arithmetic to handle values outside the expected
    range and ensure consistent wrapping.

    Parameters
    ----------
    longitude_degree_in : float
        Input longitude in degrees, expected to be in the range [-180°, 180°],
        but values outside this range are automatically wrapped. Can be any
        real number.

    Returns
    -------
    float
        Longitude in degrees in the range [0°, 360°]. The value 0° and 360°
        represent the same meridian (Prime Meridian).

    Notes
    -----
    **Conversion Formula**:

        longitude_out = (longitude_in + 360) % 360

    This formula works by:
    1. Shifting input by +360° to handle negative values
    2. Taking modulo 360 to wrap to [0°, 360°] range

    **Examples of Conversion**:

    - 0° → 0° (Prime Meridian)
    - 90° → 90° (Eastern hemisphere)
    - -90° → 270° (Western hemisphere becomes eastern)
    - 180° → 180° (International Date Line)
    - -180° → 180° (International Date Line)

    **Handling Out-of-Range Values**:

    Values outside [-180°, 180°] are automatically wrapped:
    - 270° → 270° (already in target range)
    - -200° → 160° (wraps to eastern hemisphere)
    - 400° → 40° (wraps around once)

    **Why This Conversion Is Useful**:

    The [0°, 360°] system is particularly useful for:
    - Global datasets without discontinuity at International Date Line
    - Climate models and ocean circulation data
    - Computational domains spanning the Pacific Ocean
    - Avoiding negative longitude values in certain algorithms

    Warnings
    --------
    For longitude values at exactly -180°, this function returns 180° rather
    than 0°. The ±180° meridian and 0° meridian are on opposite sides of Earth.

    Examples
    --------
    Convert typical longitude values:

    >>> convert_180_to_360(0.0)
    0.0
    >>> convert_180_to_360(90.0)
    90.0
    >>> convert_180_to_360(-90.0)
    270.0
    >>> convert_180_to_360(180.0)
    180.0
    >>> convert_180_to_360(-180.0)
    180.0

    Handle values outside expected range:

    >>> convert_180_to_360(270.0)
    270.0
    >>> convert_180_to_360(-270.0)
    90.0
    >>> convert_180_to_360(450.0)
    90.0

    Convert coordinates for Pacific-centered map:

    >>> # Los Angeles: 118°W
    >>> lon = convert_180_to_360(-118.0)
    >>> print(f"{lon}°E")
    242.0°E

    >>> # Tokyo: 140°E
    >>> lon = convert_180_to_360(140.0)
    >>> print(f"{lon}°E")
    140.0°E

    Use with global climate data:

    >>> # Convert western hemisphere longitudes for global grid
    >>> san_francisco = convert_180_to_360(-122.4)  # 122.4°W → 237.6°E
    >>> new_york = convert_180_to_360(-74.0)  # 74°W → 286°E
    >>> london = convert_180_to_360(-0.1)  # 0.1°W → 359.9°E

    See Also
    --------
    convert_360_to_180 : Inverse conversion ([0, 360] to [-180, 180])

    References
    ----------
    .. [1] IDL Coyote longitude conversion:
           http://www.idlcoyote.com/map_tips/lonconvert.html
    """
    longitude_out = (longitude_degree_in + 360.0) % 360.0
    return longitude_out
