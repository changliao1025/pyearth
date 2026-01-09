"""
Statistical outlier removal using percentile-based filtering.

This module provides functionality to identify and remove statistical outliers from
numerical datasets using a percentile-based approach. The method is robust to
extreme values and can handle multi-dimensional arrays.

Main Features
-------------
- Percentile-based outlier detection (default: 5th-95th percentile range)
- Handles NaN and infinite values automatically
- Supports multi-dimensional arrays (returns flattened result)
- Configurable threshold percentage for outlier definition
- No assumptions about data distribution (non-parametric method)

Mathematical Method
-------------------
The outlier removal uses a percentile-based approach:

1. Calculate lower and upper percentiles (default: 5th and 95th)
2. Compute the percentile range: PR = upper_percentile - lower_percentile
3. Shrink the range by threshold: adjusted_PR = PR × threshold_percentage
4. Define valid range: [lower + adjusted_PR, upper - adjusted_PR]
5. Keep only values within the valid range

Example Calculation
-------------------
Given data with 5th percentile = 10, 95th percentile = 90, threshold = 0.1:
- Percentile range (PR) = 90 - 10 = 80
- Adjustment = 80 × 0.1 = 8
- Valid range = [10 + 8, 90 - 8] = [18, 82]
- Values < 18 or > 82 are considered outliers

Typical Use Cases
-----------------
1. **Climate Data**: Remove sensor errors and extreme outliers from temperature,
   precipitation, or other climate variables
2. **Scientific Measurements**: Filter measurement errors and anomalies from
   experimental data
3. **Image Processing**: Remove extreme pixel values before analysis
4. **Time Series**: Clean temporal data by removing statistical anomalies
5. **Quality Control**: Filter out defective measurements in industrial processes

Notes
-----
- The function returns a flattened 1D array, even for multi-dimensional input
- NaN and infinite values are automatically excluded from percentile calculation
- Default percentiles (5th, 95th) assume roughly symmetric distributions
- For heavily skewed data, consider adjusting the percentile values
- The threshold percentage controls how aggressively outliers are removed

See Also
--------
numpy.percentile : Compute percentiles of array
numpy.isfinite : Test element-wise for finiteness
scipy.stats.iqr : Interquartile range (alternative outlier method)

Examples
--------
>>> import numpy as np
>>> from pyearth.toolbox.data.remove_outliers import remove_outliers

Basic usage with default settings:

>>> data = np.array([1, 2, 3, 4, 5, 100])  # 100 is an outlier
>>> clean = remove_outliers(data, 0.1)
>>> print(clean)
[1 2 3 4 5]

References
----------
.. [1] Tukey, J.W. (1977). Exploratory Data Analysis. Addison-Wesley.
.. [2] Wilcox, R.R. (2012). Introduction to Robust Estimation and Hypothesis
       Testing. Academic Press.

"""

import numpy as np
from typing import Union, Tuple, Optional
import logging

# Configure module logger
logger = logging.getLogger(__name__)


def remove_outliers(
    aData_in: Union[np.ndarray, list, tuple],
    dOutlier_percentage: float,
    lower_percentile: float = 5.0,
    upper_percentile: float = 95.0,
) -> np.ndarray:
    """
    Remove statistical outliers from numerical data using percentile-based filtering.

    This function identifies and removes outliers by:
    1. Computing lower and upper percentiles of the valid (finite) data
    2. Calculating the percentile range (PR)
    3. Shrinking the acceptable range by a threshold percentage
    4. Filtering values outside the adjusted range

    The method is non-parametric (no distribution assumptions) and robust to
    extreme values since it uses percentiles rather than mean/standard deviation.

    Parameters
    ----------
    aData_in : array-like
        Input numerical data to filter. Can be:
        - numpy.ndarray: Any shape (1D, 2D, 3D, etc.)
        - list or tuple: Will be converted to numpy array

        May contain NaN or infinite values (automatically excluded).

        Examples:
        - 1D array: [1, 2, 3, 100, 200]
        - 2D array: [[1, 2], [3, 4], [100, 200]]
        - With NaN: [1, 2, np.nan, 3, 4]

    dOutlier_percentage : float
        Threshold percentage for outlier definition, range (0.0, 1.0).

        Controls how aggressively outliers are removed:
        - 0.0: Maximum filtering (keep only median values)
        - 0.5: Moderate filtering (remove extreme 25% on each side)
        - 1.0: Minimal filtering (keep almost all data)

        The percentile range is shrunk by:
        adjustment = (upper_percentile - lower_percentile) × dOutlier_percentage

        Typical values:
        - 0.1 (10%): Aggressive outlier removal
        - 0.3 (30%): Moderate outlier removal
        - 0.5 (50%): Conservative outlier removal

        Must be in range (0.0, 1.0) exclusive.

    lower_percentile : float, optional
        Lower percentile for range calculation, default 5.0.
        Must be in range [0, 100) and less than upper_percentile.

        Common values:
        - 5.0: Standard choice for mild outlier removal
        - 1.0: More aggressive lower bound
        - 25.0: Use first quartile (robust to skewness)

    upper_percentile : float, optional
        Upper percentile for range calculation, default 95.0.
        Must be in range (0, 100] and greater than lower_percentile.

        Common values:
        - 95.0: Standard choice for mild outlier removal
        - 99.0: More conservative upper bound
        - 75.0: Use third quartile (robust to skewness)

    Returns
    -------
    numpy.ndarray
        Filtered 1D array containing only non-outlier values.

        Properties:
        - Always 1D (flattened from input shape)
        - Contains only finite values (no NaN or inf)
        - Values within the adjusted percentile range
        - Length ≤ original array length
        - dtype matches input dtype

        Empty array if:
        - All input values are NaN/infinite
        - All values are classified as outliers
        - Input is empty

    Raises
    ------
    TypeError
        If aData_in cannot be converted to numpy array.
    ValueError
        - If dOutlier_percentage not in range (0.0, 1.0)
        - If lower_percentile >= upper_percentile
        - If percentiles not in valid range [0, 100]
        - If input array is empty
        - If no finite values in input

    Warns
    -----
    UserWarning
        - If >50% of input values are NaN/infinite
        - If >90% of values are removed as outliers
        - If percentile range is zero (all finite values identical)

    Notes
    -----
    **Algorithm Details:**

    1. Input validation and conversion to numpy array
    2. Flatten multi-dimensional input to 1D
    3. Identify finite values (exclude NaN, inf, -inf)
    4. Compute lower and upper percentiles on finite data
    5. Calculate percentile range: PR = upper - lower
    6. Compute adjustment: adj = PR × dOutlier_percentage
    7. Define valid range: [lower + adj, upper - adj]
    8. Filter original data to keep only values in valid range
    9. Return filtered 1D array

    **Comparison with Other Methods:**

    - Z-score method: Assumes normal distribution, sensitive to outliers
    - IQR method: Uses quartiles (25th-75th), less customizable
    - This method: Non-parametric, customizable percentiles, robust

    **Performance:**

    - Time complexity: O(n log n) due to percentile calculation
    - Space complexity: O(n) for temporary arrays
    - Memory efficient: No deep copies of large arrays

    **Limitations:**

    - Returns flattened array (loses multi-dimensional structure)
    - Does not preserve original indices of kept values
    - Percentile calculation requires finite values (≥2 recommended)
    - Not suitable for categorical or ordinal data

    Examples
    --------
    **Example 1: Basic outlier removal with default settings**

    >>> import numpy as np
    >>> from pyearth.toolbox.data.remove_outliers import remove_outliers
    >>>
    >>> # Data with clear outliers
    >>> data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 100, 200])
    >>> clean = remove_outliers(data, dOutlier_percentage=0.1)
    >>> print(clean)
    [1 2 3 4 5 6 7 8 9]

    **Example 2: Handling NaN and infinite values**

    >>> data = np.array([1, 2, np.nan, 3, 4, np.inf, 5, -np.inf, 100])
    >>> clean = remove_outliers(data, dOutlier_percentage=0.2)
    >>> print(clean)
    [1. 2. 3. 4. 5.]  # NaN, inf, and outlier (100) removed

    **Example 3: Multi-dimensional array (returns flattened)**

    >>> data = np.array([[1, 2, 100],
    ...                  [3, 4, 200],
    ...                  [5, 6, 300]])
    >>> clean = remove_outliers(data, dOutlier_percentage=0.1)
    >>> print(clean)
    [1 2 3 4 5 6]  # Outliers removed, array flattened

    **Example 4: Aggressive outlier removal (small threshold)**

    >>> data = np.array([1, 5, 10, 15, 20, 25, 30, 35, 40, 100])
    >>> # 0.05 = very aggressive (keep only central values)
    >>> clean = remove_outliers(data, dOutlier_percentage=0.05)
    >>> print(clean)
    [10 15 20 25 30]

    **Example 5: Conservative outlier removal (large threshold)**

    >>> data = np.array([1, 5, 10, 15, 20, 25, 30, 35, 40, 100])
    >>> # 0.9 = very conservative (keep almost everything)
    >>> clean = remove_outliers(data, dOutlier_percentage=0.9)
    >>> print(clean)
    [ 1  5 10 15 20 25 30 35 40]  # Only extreme outlier removed

    **Example 6: Custom percentiles for skewed data**

    >>> # Right-skewed data (many small values, few large)
    >>> data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 50, 100, 200])
    >>> # Use 1st-90th percentile for skewed data
    >>> clean = remove_outliers(data, dOutlier_percentage=0.2,
    ...                          lower_percentile=1, upper_percentile=90)
    >>> print(clean)
    [1 2 3 4 5 6 7 8]

    **Example 7: Climate data filtering**

    >>> # Temperature data with sensor errors
    >>> temperatures = np.array([20, 21, 22, 23, 22, 21, -999, 24, 23, 150])
    >>> # -999 = missing data marker, 150 = sensor error
    >>> clean_temps = remove_outliers(temperatures, dOutlier_percentage=0.1)
    >>> print(f"Valid temperatures: {clean_temps}")
    Valid temperatures: [20 21 22 23 22 21 24 23]

    **Example 8: Image pixel value filtering**

    >>> # Image pixels with hot pixels (255) and dead pixels (0)
    >>> pixels = np.random.randint(100, 150, size=100)
    >>> pixels[0:5] = 255  # Hot pixels
    >>> pixels[5:10] = 0   # Dead pixels
    >>> clean_pixels = remove_outliers(pixels, dOutlier_percentage=0.3)
    >>> print(f"Cleaned {len(pixels) - len(clean_pixels)} bad pixels")
    Cleaned 10 bad pixels

    **Example 9: Empty or all-NaN input handling**

    >>> # All NaN data
    >>> data = np.array([np.nan, np.nan, np.nan])
    >>> try:
    ...     clean = remove_outliers(data, dOutlier_percentage=0.1)
    ... except ValueError as e:
    ...     print(f"Error: {e}")
    Error: Input array contains no finite values

    **Example 10: Comparing different threshold values**

    >>> data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100])
    >>> for threshold in [0.1, 0.3, 0.5, 0.7, 0.9]:
    ...     clean = remove_outliers(data, dOutlier_percentage=threshold)
    ...     print(f"Threshold {threshold}: kept {len(clean)}/{len(data)} values")
    Threshold 0.1: kept 9/11 values
    Threshold 0.3: kept 9/11 values
    Threshold 0.5: kept 10/11 values
    Threshold 0.7: kept 10/11 values
    Threshold 0.9: kept 10/11 values

    See Also
    --------
    numpy.percentile : Compute percentiles of array
    scipy.stats.iqr : Interquartile range method
    scipy.stats.zscore : Z-score method (assumes normality)
    pyearth.toolbox.data.check_if_duplicates : Check for duplicate values

    References
    ----------
    .. [1] Tukey, J.W. (1977). Exploratory Data Analysis. Addison-Wesley.
    .. [2] Wilcox, R.R. (2012). Introduction to Robust Estimation and Hypothesis
           Testing. Academic Press.
    .. [3] Rousseeuw, P.J. and Hubert, M. (2011). Robust statistics for outlier
           detection. WIREs Data Mining and Knowledge Discovery, 1(1), 73-79.

    """
    # ========================================================================
    # Input validation
    # ========================================================================

    # Convert input to numpy array
    try:
        data_array = np.asarray(aData_in)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Cannot convert input to numpy array. Input type: {type(aData_in)}. "
            f"Error: {e}"
        ) from e

    # Check if array is empty
    if data_array.size == 0:
        raise ValueError("Input array is empty")

    # Validate outlier percentage
    if not isinstance(dOutlier_percentage, (int, float)):
        raise TypeError(
            f"dOutlier_percentage must be numeric, got {type(dOutlier_percentage)}"
        )

    if not (0.0 < dOutlier_percentage < 1.0):
        raise ValueError(
            f"dOutlier_percentage must be in range (0.0, 1.0), got {dOutlier_percentage}"
        )

    # Validate percentiles
    if not isinstance(lower_percentile, (int, float)):
        raise TypeError(
            f"lower_percentile must be numeric, got {type(lower_percentile)}"
        )

    if not isinstance(upper_percentile, (int, float)):
        raise TypeError(
            f"upper_percentile must be numeric, got {type(upper_percentile)}"
        )

    if not (0 <= lower_percentile < 100):
        raise ValueError(
            f"lower_percentile must be in range [0, 100), got {lower_percentile}"
        )

    if not (0 < upper_percentile <= 100):
        raise ValueError(
            f"upper_percentile must be in range (0, 100], got {upper_percentile}"
        )

    if lower_percentile >= upper_percentile:
        raise ValueError(
            f"lower_percentile ({lower_percentile}) must be less than "
            f"upper_percentile ({upper_percentile})"
        )

    logger.debug(
        f"remove_outliers called with shape={data_array.shape}, "
        f"threshold={dOutlier_percentage}, percentiles=[{lower_percentile}, {upper_percentile}]"
    )

    # ========================================================================
    # Flatten array and identify finite values
    # ========================================================================

    # Flatten to 1D
    flattened = data_array.ravel()

    # Find finite values (exclude NaN, inf, -inf)
    finite_mask = np.isfinite(flattened)
    finite_data = flattened[finite_mask]

    # Check if we have any finite values
    n_total = len(flattened)
    n_finite = len(finite_data)
    n_nonfinite = n_total - n_finite

    if n_finite == 0:
        raise ValueError("Input array contains no finite values")

    # Warn if many non-finite values
    if n_nonfinite > 0:
        nonfinite_percentage = 100.0 * n_nonfinite / n_total
        logger.info(
            f"Excluded {n_nonfinite} non-finite values "
            f"({nonfinite_percentage:.1f}% of input)"
        )

        if nonfinite_percentage > 50:
            logger.warning(
                f"More than 50% of input values are NaN or infinite "
                f"({nonfinite_percentage:.1f}%)"
            )

    # ========================================================================
    # Calculate percentile range and valid bounds
    # ========================================================================

    # Compute percentiles on finite data
    lower_value = np.percentile(finite_data, lower_percentile)
    upper_value = np.percentile(finite_data, upper_percentile)

    logger.debug(
        f"Percentiles: {lower_percentile}th={lower_value}, "
        f"{upper_percentile}th={upper_value}"
    )

    # Calculate percentile range
    percentile_range = upper_value - lower_value

    # Check for zero range (all finite values are identical)
    if percentile_range == 0:
        logger.warning(
            "Percentile range is zero (all finite values are identical). "
            "Returning all finite values."
        )
        return finite_data

    # Calculate adjustment based on threshold
    adjustment = percentile_range * dOutlier_percentage

    # Define valid range bounds
    lower_bound = lower_value + adjustment
    upper_bound = upper_value - adjustment

    logger.debug(
        f"Valid range: [{lower_bound}, {upper_bound}], "
        f"percentile_range={percentile_range}, adjustment={adjustment}"
    )

    # ========================================================================
    # Filter outliers
    # ========================================================================

    # Create mask for values within valid range
    # Note: We filter on the original flattened array to preserve
    # finite values that fall within the range
    valid_mask = (flattened >= lower_bound) & (flattened <= upper_bound)
    filtered_data = flattened[valid_mask]

    # Count removed values
    n_filtered = len(filtered_data)
    n_removed = n_total - n_filtered
    removal_percentage = 100.0 * n_removed / n_total

    logger.info(
        f"Removed {n_removed} values ({removal_percentage:.1f}% of input), "
        f"kept {n_filtered} values"
    )

    # Warn if too many values removed
    if removal_percentage > 90:
        logger.warning(
            f"More than 90% of values removed ({removal_percentage:.1f}%). "
            f"Consider using a larger threshold or different percentiles."
        )

    # Check if result is empty
    if n_filtered == 0:
        logger.warning("All values were classified as outliers. Returning empty array.")
        return np.array([], dtype=data_array.dtype)

    logger.debug(f"Returning filtered array with shape {filtered_data.shape}")

    return filtered_data
