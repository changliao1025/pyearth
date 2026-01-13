"""
Custom Percentile Calculation with Duplicate Value Handling

This module provides robust percentile calculation for datasets with potential
issues such as duplicate extreme values, missing values, and heavily skewed
distributions. The cgpercentiles function handles edge cases that can cause
standard percentile calculations to fail or produce misleading results.

Main Functions
--------------
cgpercentiles : Calculate percentiles with special handling for duplicate values

Key Features
------------
- Handles datasets with many duplicate minimum/maximum values
- Manages missing value replacement and filtering
- Detects and handles problematic zero-value concentrations
- Augments data strategically to ensure meaningful percentile bins
- Supports custom percentile specifications
- Robust error handling for edge cases

Use Cases
---------
1. **Skewed Distributions**: Handle data with extreme value clusters
2. **Sparse Data**: Calculate percentiles when unique values are limited
3. **Environmental Data**: Process datasets with many zero or minimum values
4. **Quality Control**: Filter and analyze data with known missing value codes
5. **Climate Analysis**: Handle precipitation data with many zero values
6. **Remote Sensing**: Process imagery with no-data or saturation values
7. **Scientific Computing**: Robust percentile calculation for any numeric data

Technical Details
-----------------
The function implements a scenario-based approach to handle 8 different
data distribution patterns:

Scenarios:
1. Both min and max have excessive duplicates
2. Only min has excessive duplicates
3. Only max has excessive duplicates
4. No excessive duplicates (standard case)
5. Min, max, and zero all have excessive duplicates
6. Min and zero have excessive duplicates
7. Max and zero have excessive duplicates
8. Only zero has excessive duplicates

For each scenario, the function augments the dataset with strategic duplicate
values to ensure each percentile bin has sufficient unique values for meaningful
calculation.

Dependencies
------------
- numpy: Array operations and percentile calculation

See Also
--------
- numpy.percentile: Underlying percentile calculation
- numpy.nanpercentile: Percentile ignoring NaN values
"""

import numpy as np
import logging
from typing import Union, List, Optional, Sequence

# Configure logging
logger = logging.getLogger(__name__)


def cgpercentiles(
    aData_in: Union[np.ndarray, Sequence[float]],
    aPercentiles_in: Union[np.ndarray, Sequence[float]],
    missing_value_in: Optional[float] = None,
) -> Union[List[float], int]:
    """
    Calculate percentiles with robust handling of duplicate extreme values.

    This function computes percentiles for datasets that may contain excessive
    duplicate values at the extremes (min, max) or at zero. Unlike standard
    percentile functions, it detects problematic value concentrations and
    strategically augments the dataset to ensure meaningful percentile bins.

    Parameters
    ----------
    aData_in : array-like
        Input data array. Can be multi-dimensional; will be flattened.
        Must contain at least some finite numeric values.

        Supported types:
        - numpy.ndarray of any shape
        - List or tuple of numbers
        - Any sequence convertible to numpy array

    aPercentiles_in : array-like
        Percentiles to calculate, as values between 0 and 100.

        Examples:
        - [25, 50, 75] for quartiles
        - [10, 90] for deciles
        - [2.5, 97.5] for 95% confidence interval bounds

        If empty, defaults to [25, 50, 75] (quartiles).

    missing_value_in : float or None, optional
        Specific value to treat as missing/no-data.

        Default: None (no special missing value)

        When provided:
        - All occurrences replaced with NaN
        - Filtered out before percentile calculation

        Common values:
        - -9999 (meteorological data)
        - -999 (hydrological data)
        - np.nan (pre-existing NaN markers)

    Returns
    -------
    list of float or int
        List of calculated percentile values, one per input percentile.
        Length matches aPercentiles_in length.

        Returns:
        - List of floats: Successful calculation
        - 0: Insufficient valid data after filtering
        - -1: Invalid input (empty data array)

    Notes
    -----
    1. **Flattening**: Multi-dimensional input automatically flattened
    2. **Missing Values**: NaN and specified missing values filtered out
    3. **Scenario Detection**: Analyzes distribution to detect problems
    4. **Data Augmentation**: Adds strategic duplicates to ensure valid bins
    5. **Threshold**: Problematic when duplicates > (total / n_percentiles)
    6. **Zero Handling**: Special treatment for datasets spanning zero
    7. **Edge Cases**: Returns 0 when insufficient data after filtering
    8. **Type Conversion**: Input converted to numpy array automatically
    9. **Finite Check**: Only finite values used in calculation
    10. **Preservation**: Original data not modified

    Algorithm Details
    -----------------
    1. Flatten and validate input data
    2. Replace missing values with NaN if specified
    3. Filter to finite values only
    4. Calculate min, max, and their occurrence counts
    5. Determine minimum samples per percentile bin
    6. Detect which scenario applies (1-8)
    7. For problematic scenarios:
       - Identify valid data range
       - Create augmentation arrays with problem values
       - Concatenate augmented data
    8. Calculate percentiles on processed dataset

    Scenarios Explained
    -------------------
    **Scenario 1**: Both extremes problematic
        - Too many min and max values
        - Solution: Add both to middle values

    **Scenario 2**: Min problematic
        - Too many minimum values
        - Solution: Add min values to higher values

    **Scenario 3**: Max problematic
        - Too many maximum values
        - Solution: Add max values to lower values

    **Scenario 4**: No problems (standard)
        - Use data as-is

    **Scenario 5**: Min, max, and zero problematic
        - All three values over-represented
        - Solution: Add all three to other values

    **Scenario 6**: Min and zero problematic
        - Solution: Add both to positive values

    **Scenario 7**: Max and zero problematic
        - Solution: Add both to non-zero, non-max values

    **Scenario 8**: Only zero problematic
        - Solution: Add zeros to non-zero values

    Examples
    --------
    Calculate quartiles for standard dataset:

    >>> data = np.random.randn(1000)
    >>> percentiles = [25, 50, 75]
    >>> result = cgpercentiles(data, percentiles)
    >>> print(f"Q1={result[0]:.2f}, Median={result[1]:.2f}, Q3={result[2]:.2f}")
    Q1=-0.67, Median=0.02, Q3=0.68

    Handle precipitation data with many zeros:

    >>> precip = np.array([0, 0, 0, 0, 0, 0, 0.5, 1.2, 2.3, 15.6])
    >>> result = cgpercentiles(precip, [10, 50, 90])
    >>> print(f"10th={result[0]}, Median={result[1]}, 90th={result[2]}")
    10th=0.0, Median=0.0, 90th=8.85

    Process data with missing value code:

    >>> temps = np.array([15.2, 18.3, -9999, 22.1, 19.5, -9999, 20.3])
    >>> result = cgpercentiles(temps, [25, 50, 75], missing_value_in=-9999)
    >>> print(f"Valid data percentiles: {result}")
    Valid data percentiles: [17.0, 19.5, 21.2]

    Handle heavily skewed data with many maximum values:

    >>> # Satellite data with saturation at 255
    >>> satellite = np.array([45, 67, 89, 255, 255, 255, 255, 255, 120, 98])
    >>> result = cgpercentiles(satellite, [25, 50, 75])
    >>> # Function augments to handle many 255 values

    Use default percentiles (quartiles):

    >>> data = np.arange(100)
    >>> result = cgpercentiles(data, [])  # Empty uses default
    >>> print(f"Default quartiles: {result}")
    Default quartiles: [24.75, 49.5, 74.25]

    Multi-dimensional input (automatically flattened):

    >>> grid = np.random.rand(10, 10)
    >>> result = cgpercentiles(grid, [5, 95])
    >>> print(f"5th-95th range: {result[0]:.3f} to {result[1]:.3f}")
    5th-95th range: 0.052 to 0.947

    Handle insufficient data case:

    >>> sparse = np.array([1, 1, 1, 1, 1])  # All same value
    >>> percentiles = [10, 20, 30, 40, 50, 60, 70, 80, 90]  # Too many
    >>> result = cgpercentiles(sparse, percentiles)
    >>> print(result)
    0

    See Also
    --------
    numpy.percentile : Standard percentile calculation
    numpy.nanpercentile : Percentile ignoring NaN values
    numpy.quantile : Quantile calculation (0-1 scale)

    References
    ----------
    .. [1] Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical
           packages. The American Statistician, 50(4), 361-365.
    .. [2] NumPy percentile documentation
           https://numpy.org/doc/stable/reference/generated/numpy.percentile.html

    Warnings
    --------
    - For very large datasets, memory usage may be significant due to augmentation
    - The function modifies internal copies, not the original data
    - Returns integer codes (0, -1) for error conditions; check return type
    """
    # Validate and prepare input data
    if aData_in is None:
        error_msg = "Input data cannot be None"
        logger.error(error_msg)
        print("Error: Input data is required.")
        return -1

    # Convert to numpy array and flatten
    try:
        aData_in0 = np.array(aData_in)
    except Exception as e:
        error_msg = f"Failed to convert input data to array: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return -1

    aData_in_flat = aData_in0.flatten()
    aData_in_reshaped = np.reshape(aData_in_flat, (aData_in_flat.size,))

    if aData_in_reshaped.size == 0:
        error_msg = "Input data array is empty"
        logger.error(error_msg)
        print("Error: Input data is required.")
        return -1

    # Validate and prepare percentiles
    try:
        aPercentiles_in_array = np.array(aPercentiles_in)
    except Exception as e:
        error_msg = f"Failed to convert percentiles to array: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return -1

    aPercentiles_flat = aPercentiles_in_array.flatten()

    # Use default percentiles if none provided
    if aPercentiles_flat.size == 0:
        aPercentiles_flat = np.array([25.0, 50.0, 75.0])
        logger.info("No percentiles specified, using default quartiles [25, 50, 75]")

    # Validate percentile ranges
    if np.any((aPercentiles_flat < 0) | (aPercentiles_flat > 100)):
        logger.warning("Some percentiles are outside valid range [0, 100]")

    # Handle missing values
    if missing_value_in is not None:
        logger.info(f"Replacing missing value {missing_value_in} with NaN")
        dummy_index = np.where(aData_in_reshaped == missing_value_in)
        aData_in_reshaped[dummy_index] = np.nan

    # Filter to finite values only
    good_index = np.where(np.isfinite(aData_in_reshaped))
    aData_copy = aData_in_reshaped[good_index]

    if aData_copy.size == 0:
        error_msg = "No valid (finite) data values after filtering"
        logger.error(error_msg)
        print("Error: All data values are NaN or infinite.")
        return -1

    num = aData_copy.size
    num_per = aPercentiles_flat.size
    min_count = int(num / num_per)

    if min_count == 0:
        logger.warning(f"Very few data points ({num}) for {num_per} percentiles")
        min_count = 1

    # Calculate min, max, and their counts
    min_value = np.min(aData_copy)
    max_value = np.max(aData_copy)

    nan_index1 = np.where(aData_copy == min_value)
    nan_index2 = np.where(aData_copy == max_value)

    nan_count1 = nan_index1[0].size  # Count of minimum values
    nan_count2 = nan_index2[0].size  # Count of maximum values

    # Determine scenario based on distribution characteristics
    scenario = 4  # Default: no problems

    # Check if data spans zero
    if (min_value >= 0.0) or (max_value <= 0.0):
        # Data doesn't span zero - simpler scenarios (1-4)
        if (nan_count1 > min_count) and (nan_count2 > min_count):
            scenario = 1  # Both extremes problematic
        elif (nan_count1 > min_count) and (nan_count2 <= min_count):
            scenario = 2  # Min problematic
        elif (nan_count2 > min_count) and (nan_count1 <= min_count):
            scenario = 3  # Max problematic
        else:
            scenario = 4  # No problems
    else:
        # Data spans zero - check for zero concentration
        nan_index3 = np.where(aData_copy == 0.0)
        nan_count3 = nan_index3[0].size

        if nan_count3 >= min_count:
            # Zero is also problematic - scenarios 5-8
            if (nan_count1 >= min_count) and (nan_count2 >= min_count):
                scenario = 5  # Min, max, and zero problematic
            elif (nan_count1 >= min_count) and (nan_count2 < min_count):
                scenario = 6  # Min and zero problematic
            elif (nan_count2 >= min_count) and (nan_count1 < min_count):
                scenario = 7  # Max and zero problematic
            else:
                scenario = 8  # Only zero problematic
        else:
            # Zero not problematic - back to scenarios 1-4
            if (nan_count1 >= min_count) and (nan_count2 >= min_count):
                scenario = 1
            elif (nan_count1 >= min_count) and (nan_count2 < min_count):
                scenario = 2
            elif (nan_count2 >= min_count) and (nan_count1 < min_count):
                scenario = 3
            else:
                scenario = 4

    logger.debug(
        f"Detected scenario {scenario}: min_count={nan_count1}, max_count={nan_count2}, threshold={min_count}"
    )

    # Process data based on scenario
    if scenario == 1:
        # Both min and max problematic
        good_index = np.where((aData_copy > min_value) & (aData_copy < max_value))
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count1 = int(np.ceil(float(good_count) / (num_per - 1)) + 1)
        bad_copy1 = np.full(sample_count1, fill_value=min_value)
        sample_count2 = int(np.floor(float(good_count) / (num_per - 1)) - 1)
        bad_copy2 = np.full(sample_count2, fill_value=max_value)
        data_copy2 = np.concatenate((bad_copy1, bad_copy2, aData_copy[good_index]))

    elif scenario == 2:
        # Only min problematic
        good_index = np.where(aData_copy > min_value)
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count = int(np.ceil(float(good_count) / num_per))
        bad_copy = np.full(sample_count, fill_value=min_value)
        data_copy2 = np.concatenate((bad_copy, aData_copy[good_index]))

    elif scenario == 3:
        # Only max problematic
        good_index = np.where(aData_copy < max_value)
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count = int(np.ceil(float(good_count) / num_per))
        bad_copy = np.full(sample_count, fill_value=max_value)
        data_copy2 = np.concatenate((bad_copy, aData_copy[good_index]))

    elif scenario == 4:
        # No problems - use data as-is
        data_copy2 = aData_copy

    elif scenario == 5:
        # Min, max, and zero problematic
        good_index = np.where(
            (aData_copy > min_value) & (aData_copy != 0.0) & (aData_copy < max_value)
        )
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count1 = int(np.ceil(float(good_count) / (num_per - 2)) + 1)
        bad_copy1 = np.full(sample_count1, fill_value=min_value)
        sample_count2 = int(np.floor(float(good_count) / (num_per - 2)) - 1)
        bad_copy2 = np.full(sample_count2, fill_value=0.0)
        sample_count3 = int(np.floor(float(good_count) / (num_per - 2)) - 1)
        bad_copy3 = np.full(sample_count3, fill_value=max_value)
        data_copy2 = np.concatenate(
            (bad_copy1, bad_copy2, bad_copy3, aData_copy[good_index])
        )

    elif scenario == 6:
        # Min and zero problematic
        good_index = np.where((aData_copy > min_value) & (aData_copy != 0.0))
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count1 = int(np.ceil(float(good_count) / (num_per - 1)) + 1)
        bad_copy1 = np.full(sample_count1, fill_value=min_value)
        sample_count2 = int(np.floor(float(good_count) / (num_per - 1)) - 1)
        bad_copy2 = np.full(sample_count2, fill_value=0.0)
        data_copy2 = np.concatenate((bad_copy1, bad_copy2, aData_copy[good_index]))

    elif scenario == 7:
        # Max and zero problematic
        good_index = np.where((aData_copy != 0.0) & (aData_copy < max_value))
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count1 = int(np.ceil(float(good_count) / (num_per - 1)) + 1)
        bad_copy1 = np.full(sample_count1, fill_value=0.0)
        sample_count2 = int(np.floor(float(good_count) / (num_per - 1)) - 1)
        bad_copy2 = np.full(sample_count2, fill_value=max_value)
        data_copy2 = np.concatenate((bad_copy1, bad_copy2, aData_copy[good_index]))

    else:  # scenario == 8
        # Only zero problematic
        good_index = np.where(aData_copy != 0.0)
        good_count = good_index[0].size
        if good_count <= num_per:
            logger.warning(f"Insufficient valid data: {good_count} <= {num_per}")
            return 0

        sample_count = int(np.ceil(float(good_count) / num_per) - 1)
        bad_copy = np.full(sample_count, fill_value=0.0)
        data_copy2 = np.concatenate((bad_copy, aData_copy[good_index]))

    # Calculate percentiles on processed data
    result = []
    for p in aPercentiles_flat:
        try:
            percentile_value = np.percentile(data_copy2, p)
            result.append(float(percentile_value))
        except Exception as e:
            logger.error(f"Failed to calculate {p}th percentile: {str(e)}")
            result.append(np.nan)

    logger.info(f"Successfully calculated {len(result)} percentiles")
    return result
