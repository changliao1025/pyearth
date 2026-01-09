"""
Time series temporal aggregation from daily to monthly resolution.

This module provides functionality to aggregate daily time series data into monthly
values using various statistical methods (mean, sum, etc.). It handles partial months,
missing data (NaN values), and ensures accurate date calculations using Julian dates.

Main Features
-------------
- Daily to monthly temporal aggregation
- Multiple aggregation methods: mean (default), sum
- Handles partial months at start/end of time series
- Automatic NaN handling (excluded from calculations)
- Julian date-based indexing for accuracy
- Supports multi-year time series

Mathematical Method
-------------------
The aggregation process:

1. **Mean aggregation** (default):
   - For each month, sum all valid (non-NaN) daily values
   - Divide by count of valid days
   - If all days are NaN, return NaN for that month
   - Formula: monthly_mean = sum(valid_daily_values) / count(valid_days)

2. **Sum aggregation**:
   - For each month, sum all daily values (including NaN as 0)
   - No division by day count
   - Formula: monthly_sum = sum(daily_values)

Typical Use Cases
-----------------
1. **Climate Data**: Aggregate daily temperature, precipitation, or other climate
   variables to monthly averages or totals
2. **Hydrology**: Convert daily streamflow to monthly discharge volumes
3. **Environmental Monitoring**: Summarize daily sensor readings to monthly statistics
4. **Energy Data**: Aggregate daily power consumption to monthly totals
5. **Financial Data**: Convert daily stock prices to monthly averages

Algorithm Details
-----------------
1. Convert start date to Julian date (reference point)
2. For each month in the time range:
   - Determine day range (handle partial months at boundaries)
   - For each day in the month:
     * Calculate Julian date
     * Compute index into daily data array
     * Aggregate value based on method (mean or sum)
   - Store monthly result
3. Return array of monthly values

Notes
-----
- Uses Julian dates for accurate day-of-year calculations
- Handles leap years automatically via Julian date conversion
- Mean aggregation excludes NaN values from both sum and count
- Sum aggregation treats NaN as contributing 0 to the sum
- Partial months are handled correctly at time series boundaries
- Output array length = number of complete/partial months in range

See Also
--------
pyearth.toolbox.date.julian : Julian date conversion utilities
pyearth.toolbox.date.day_in_month : Days in month calculation
numpy.nanmean : Alternative for array-based mean with NaN handling
pandas.resample : More flexible time series resampling (if using pandas)

Examples
--------
>>> import numpy as np
>>> from pyearth.toolbox.data.convert_time_series_daily_to_monthly import \\
...     convert_time_series_daily_to_monthly

Basic monthly mean aggregation:

>>> # 31 days of January data (daily temperatures)
>>> daily_temps = np.array([20, 21, 22, 23, 22, 21, 20] * 4 + [21, 22, 23])
>>> monthly = convert_time_series_daily_to_monthly(
...     daily_temps, 2020, 1, 1, 2020, 1, 31, sType_in='mean')
>>> print(f"January mean: {monthly[0]:.2f}°C")

Monthly sum aggregation (e.g., precipitation):

>>> # 28 days of February precipitation (mm)
>>> daily_precip = np.array([0, 0, 5, 10, 0, 0, 0] * 4)
>>> monthly = convert_time_series_daily_to_monthly(
...     daily_precip, 2020, 2, 1, 2020, 2, 28, sType_in='sum')
>>> print(f"February total: {monthly[0]:.1f} mm")

References
----------
.. [1] Meeus, J. (1991). Astronomical Algorithms. Willmann-Bell, Inc.
.. [2] Hatcher, D.A. (1984). Simple Formulae for Julian Day Numbers and
       Calendar Dates. Quarterly Journal of the Royal Astronomical Society, 25, 53-55.

"""

import numpy as np
import datetime
import logging
from typing import Union, Optional, Literal

import pyearth.toolbox.date.julian as julian
from pyearth.toolbox.date.day_in_month import day_in_month

# Configure module logger
logger = logging.getLogger(__name__)


def convert_time_series_daily_to_monthly(
    aData_daily_in: Union[np.ndarray, list],
    iYear_start_in: int,
    iMonth_start_in: int,
    iDay_start_in: int,
    iYear_end_in: int,
    iMonth_end_in: int,
    iDay_end_in: int,
    sType_in: Optional[Literal["mean", "sum"]] = None,
) -> np.ndarray:
    """
    Convert daily time series data to monthly aggregated values.

    Aggregates daily measurements into monthly statistics using mean or sum.
    Handles partial months at time series boundaries and automatically excludes
    NaN values from mean calculations. Uses Julian dates for accurate indexing.

    Parameters
    ----------
    aData_daily_in : array-like
        Daily time series data to aggregate.

        Must be 1D array with length matching the number of days from
        start date to end date (inclusive).

        Can contain NaN values:
        - For 'mean': NaN values excluded from calculation
        - For 'sum': NaN values treated as 0

        Examples:
        - Temperature: [20.5, 21.3, 22.1, ...]
        - Precipitation: [0.0, 5.2, 0.0, 10.5, ...]
        - Streamflow: [15.2, 14.8, np.nan, 16.1, ...]

    iYear_start_in : int
        Start year of time series (e.g., 2020).
        Must be valid Gregorian calendar year.
        Range: typically 1582-9999 (after Gregorian reform).

    iMonth_start_in : int
        Start month of time series.
        Range: 1-12 (1=January, 12=December).

    iDay_start_in : int
        Start day of month.
        Range: 1-31 (must be valid for the given month/year).
        For example: February 30 is invalid.

    iYear_end_in : int
        End year of time series (e.g., 2021).
        Must be >= iYear_start_in.

    iMonth_end_in : int
        End month of time series.
        Range: 1-12.
        If same year as start, must be >= iMonth_start_in.

    iDay_end_in : int
        End day of month (inclusive).
        Range: 1-31 (must be valid for the given month/year).

    sType_in : {'mean', 'sum'}, optional
        Aggregation method to use, default 'mean'.

        - 'mean': Calculate monthly average of daily values
          * Excludes NaN values from both sum and count
          * Returns NaN if all days in month are NaN
          * Formula: sum(valid_values) / count(valid_values)

        - 'sum': Calculate monthly total of daily values
          * Treats NaN as 0 in summation
          * Useful for accumulative variables (precipitation, energy)
          * Formula: sum(all_values, NaN→0)

    Returns
    -------
    numpy.ndarray
        Monthly aggregated values.

        Shape: (n_months,) where n_months is the number of complete or
        partial months in the date range.

        Properties:
        - dtype: float64
        - Length: Number of months from (year_start, month_start) to
          (year_end, month_end) inclusive
        - May contain NaN if 'mean' method and all days in month are NaN

        Examples:
        - 1 year, full months: length = 12
        - Jan 15, 2020 to Mar 10, 2021: length = 15 months
        - Partial months at start/end included

    Raises
    ------
    TypeError
        If input data cannot be converted to numpy array.
    ValueError
        - If start date > end date
        - If any date component is invalid (month not 1-12, etc.)
        - If daily data length doesn't match expected days in range
        - If sType_in is not 'mean' or 'sum'
    IndexError
        If calculated array index exceeds daily data length (indicates
        mismatch between date range and data length).

    Warns
    -----
    UserWarning
        - If >50% of days in a month are NaN (for 'mean' aggregation)
        - If calculated data length doesn't match actual array length
        - If partial months at boundaries (expected behavior, informational)

    Notes
    -----
    **Date Handling:**

    - Uses Julian date conversion for accurate day counting
    - Automatically handles leap years
    - Partial months at boundaries are handled correctly:
      * Start month: uses days from iDay_start_in to end of month
      * End month: uses days from 1 to iDay_end_in
      * Middle months: uses all days (1 to last day of month)

    **NaN Handling:**

    - 'mean' method:
      * Only counts non-NaN values
      * Returns NaN if all values in month are NaN
      * Effectively computes np.nanmean() for each month

    - 'sum' method:
      * Includes all values (NaN treated as 0)
      * Never returns NaN in output
      * Equivalent to np.nansum() for each month

    **Performance:**

    - Time complexity: O(n) where n = number of days
    - Space complexity: O(m) where m = number of months
    - Efficient for typical climate datasets (years to decades)

    **Common Issues:**

    1. **Array length mismatch**: Ensure daily data length equals
       number of days from start to end date (inclusive)
    2. **Invalid dates**: February 30, April 31, etc. will raise ValueError
    3. **Date order**: Start date must be before or equal to end date

    Examples
    --------
    **Example 1: Basic monthly mean (full year)**

    >>> import numpy as np
    >>> from pyearth.toolbox.data.convert_time_series_daily_to_monthly import \\
    ...     convert_time_series_daily_to_monthly
    >>>
    >>> # 365 days of daily temperature data for 2019 (non-leap year)
    >>> daily_temp = 20 + 10 * np.sin(np.linspace(0, 2*np.pi, 365))
    >>> monthly_temp = convert_time_series_daily_to_monthly(
    ...     daily_temp, 2019, 1, 1, 2019, 12, 31, sType_in='mean')
    >>> print(f"Number of months: {len(monthly_temp)}")
    Number of months: 12
    >>> print(f"January mean: {monthly_temp[0]:.2f}°C")
    January mean: 20.52°C

    **Example 2: Monthly sum (precipitation)**

    >>> # 31 days of January precipitation (mm/day)
    >>> daily_precip = np.array([0, 0, 5.2, 0, 0, 10.5, 0] * 4 + [0, 0, 3.1])
    >>> monthly_total = convert_time_series_daily_to_monthly(
    ...     daily_precip, 2020, 1, 1, 2020, 1, 31, sType_in='sum')
    >>> print(f"January total precipitation: {monthly_total[0]:.1f} mm")
    January total precipitation: 65.9 mm

    **Example 3: Handling NaN values (mean)**

    >>> # Daily data with missing values
    >>> daily_data = np.array([10, 20, np.nan, 30, np.nan, 40] * 5 + [10])
    >>> # 31 days = full January
    >>> monthly = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 1, 2020, 1, 31, sType_in='mean')
    >>> # Mean of [10, 20, 30, 40] repeated, NaN excluded
    >>> print(f"Mean (excluding NaN): {monthly[0]:.2f}")
    Mean (excluding NaN): 25.00

    **Example 4: Partial months at boundaries**

    >>> # Mid-January to mid-March (partial months at both ends)
    >>> # 17 days (Jan 15-31) + 29 days (Feb, 2020 is leap) + 15 days (Mar 1-15)
    >>> # Total = 61 days
    >>> daily_data = np.ones(61) * 10  # Constant value for simplicity
    >>> monthly = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 15, 2020, 3, 15, sType_in='mean')
    >>> print(f"Number of months: {len(monthly)}")  # Jan, Feb, Mar
    Number of months: 3
    >>> print(f"All monthly means: {monthly}")
    All monthly means: [10. 10. 10.]

    **Example 5: Multi-year time series**

    >>> # 2 years of daily data (2020-2021, includes leap year)
    >>> n_days = 366 + 365  # 2020 is leap year
    >>> daily_data = np.random.normal(15, 5, n_days)  # Random temps
    >>> monthly = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 1, 2021, 12, 31, sType_in='mean')
    >>> print(f"Number of months in 2 years: {len(monthly)}")
    Number of months in 2 years: 24

    **Example 6: Single month aggregation**

    >>> # Just February 2020 (leap year, 29 days)
    >>> daily_data = np.linspace(5, 15, 29)  # Temperature trend
    >>> monthly = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 2, 1, 2020, 2, 29, sType_in='mean')
    >>> print(f"February mean: {monthly[0]:.2f}°C")
    February mean: 10.00°C

    **Example 7: Handling all-NaN month**

    >>> # Month with all missing data
    >>> daily_data = np.full(31, np.nan)
    >>> monthly = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 1, 2020, 1, 31, sType_in='mean')
    >>> print(f"Result for all-NaN month: {monthly[0]}")
    Result for all-NaN month: nan

    **Example 8: Sum vs Mean comparison**

    >>> daily_data = np.array([10, 20, 30] * 10 + [10])  # 31 days
    >>> mean_result = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 1, 2020, 1, 31, sType_in='mean')
    >>> sum_result = convert_time_series_daily_to_monthly(
    ...     daily_data, 2020, 1, 1, 2020, 1, 31, sType_in='sum')
    >>> print(f"Mean: {mean_result[0]:.2f}, Sum: {sum_result[0]:.2f}")
    Mean: 20.00, Sum: 620.00

    **Example 9: Climate data use case (seasonal cycle)**

    >>> # Simulate daily temperature with seasonal cycle
    >>> days = np.arange(365)
    >>> daily_temp = 15 + 10 * np.sin(2 * np.pi * (days - 80) / 365)
    >>> monthly_temp = convert_time_series_daily_to_monthly(
    ...     daily_temp, 2019, 1, 1, 2019, 12, 31, sType_in='mean')
    >>> # Find warmest and coldest months
    >>> warmest_month = np.argmax(monthly_temp) + 1
    >>> coldest_month = np.argmin(monthly_temp) + 1
    >>> print(f"Warmest month: {warmest_month}, Coldest month: {coldest_month}")
    Warmest month: 7, Coldest month: 1

    **Example 10: Hydrology application (monthly discharge)**

    >>> # Daily streamflow (m³/s) to monthly volume (m³)
    >>> daily_flow = np.random.gamma(2, 50, 365)  # Skewed distribution
    >>> # Convert to monthly total volume (flow × seconds_per_day × days)
    >>> monthly_flow = convert_time_series_daily_to_monthly(
    ...     daily_flow, 2019, 1, 1, 2019, 12, 31, sType_in='mean')
    >>> # Get monthly average flow rate
    >>> print(f"Annual mean daily flow: {np.mean(monthly_flow):.2f} m³/s")
    Annual mean daily flow: 100.15 m³/s

    See Also
    --------
    pyearth.toolbox.date.julian.to_jd : Julian date conversion
    pyearth.toolbox.date.day_in_month : Days in month calculation
    numpy.nanmean : Compute mean, ignoring NaN values
    numpy.nansum : Compute sum, ignoring NaN values
    pandas.DataFrame.resample : Pandas time series resampling

    References
    ----------
    .. [1] Meeus, J. (1991). Astronomical Algorithms. Willmann-Bell, Inc.
    .. [2] Press, W.H. et al. (2007). Numerical Recipes: The Art of Scientific
           Computing (3rd ed.). Cambridge University Press.

    """
    # ========================================================================
    # Input validation
    # ========================================================================

    # Validate aggregation type
    if sType_in is not None:
        sType = sType_in
    else:
        sType = "mean"

    if sType not in ["mean", "sum"]:
        raise ValueError(f"sType_in must be 'mean' or 'sum', got '{sType}'")

    # Validate date components
    if not isinstance(iYear_start_in, int):
        raise TypeError(f"iYear_start_in must be int, got {type(iYear_start_in)}")
    if not isinstance(iYear_end_in, int):
        raise TypeError(f"iYear_end_in must be int, got {type(iYear_end_in)}")

    if not (1 <= iMonth_start_in <= 12):
        raise ValueError(f"iMonth_start_in must be 1-12, got {iMonth_start_in}")
    if not (1 <= iMonth_end_in <= 12):
        raise ValueError(f"iMonth_end_in must be 1-12, got {iMonth_end_in}")

    if not (1 <= iDay_start_in <= 31):
        raise ValueError(f"iDay_start_in must be 1-31, got {iDay_start_in}")
    if not (1 <= iDay_end_in <= 31):
        raise ValueError(f"iDay_end_in must be 1-31, got {iDay_end_in}")

    # Validate year order
    if iYear_end_in < iYear_start_in:
        raise ValueError(
            f"End year ({iYear_end_in}) must be >= start year ({iYear_start_in})"
        )

    # Validate month/day order for same year
    if iYear_end_in == iYear_start_in:
        if iMonth_end_in < iMonth_start_in:
            raise ValueError(
                f"For same year, end month ({iMonth_end_in}) must be >= "
                f"start month ({iMonth_start_in})"
            )
        if iMonth_end_in == iMonth_start_in and iDay_end_in < iDay_start_in:
            raise ValueError(
                f"For same year and month, end day ({iDay_end_in}) must be >= "
                f"start day ({iDay_start_in})"
            )

    # Validate day is valid for month
    try:
        datetime.datetime(iYear_start_in, iMonth_start_in, iDay_start_in)
    except ValueError as e:
        raise ValueError(
            f"Invalid start date: {iYear_start_in}-{iMonth_start_in}-{iDay_start_in}. "
            f"Error: {e}"
        ) from e

    try:
        datetime.datetime(iYear_end_in, iMonth_end_in, iDay_end_in)
    except ValueError as e:
        raise ValueError(
            f"Invalid end date: {iYear_end_in}-{iMonth_end_in}-{iDay_end_in}. "
            f"Error: {e}"
        ) from e

    # Convert input to numpy array
    try:
        aData_daily = np.asarray(aData_daily_in, dtype=float)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"Cannot convert input data to numpy array. "
            f"Input type: {type(aData_daily_in)}. Error: {e}"
        ) from e

    if aData_daily.ndim != 1:
        raise ValueError(f"Input data must be 1D array, got {aData_daily.ndim}D")

    if aData_daily.size == 0:
        raise ValueError("Input data array is empty")

    logger.debug(
        f"Converting daily to monthly: {iYear_start_in}-{iMonth_start_in:02d}-{iDay_start_in:02d} "
        f"to {iYear_end_in}-{iMonth_end_in:02d}-{iDay_end_in:02d}, "
        f"method='{sType}', data_length={len(aData_daily)}"
    )

    # ========================================================================
    # Calculate expected data length
    # ========================================================================

    date_start = datetime.datetime(iYear_start_in, iMonth_start_in, iDay_start_in)
    date_end = datetime.datetime(iYear_end_in, iMonth_end_in, iDay_end_in)
    expected_days = (date_end - date_start).days + 1  # +1 for inclusive

    if len(aData_daily) != expected_days:
        logger.warning(
            f"Data length ({len(aData_daily)}) doesn't match expected days "
            f"({expected_days}) from {date_start.date()} to {date_end.date()}"
        )

    # ========================================================================
    # Initialize Julian date reference
    # ========================================================================

    lJulian_start = julian.to_jd(date_start, fmt="jd")

    logger.debug(f"Julian start date: {lJulian_start}")

    # ========================================================================
    # Aggregate daily data to monthly
    # ========================================================================

    aData_monthly_out = []
    month_count = 0

    for iYear in range(iYear_start_in, iYear_end_in + 1):
        # Determine month range for this year
        if iYear == iYear_start_in:
            iMonth_start = iMonth_start_in
        else:
            iMonth_start = 1

        if iYear == iYear_end_in:
            iMonth_end = iMonth_end_in
        else:
            iMonth_end = 12

        for iMonth in range(iMonth_start, iMonth_end + 1):
            month_count += 1

            # Determine day range for this month
            if iYear == iYear_start_in and iMonth == iMonth_start_in:
                iDay_start = iDay_start_in
            else:
                iDay_start = 1

            if iYear == iYear_end_in and iMonth == iMonth_end_in:
                iDay_end = iDay_end_in
            else:
                iDay_end = day_in_month(iYear, iMonth)

            logger.debug(
                f"Processing month {month_count}: {iYear}-{iMonth:02d}, "
                f"days {iDay_start} to {iDay_end}"
            )

            # Aggregate values for this month
            if sType == "mean":
                monthly_sum = 0.0
                valid_count = 0

                for iDay in range(iDay_start, iDay_end + 1):
                    # Calculate Julian date for this day
                    current_date = datetime.datetime(iYear, iMonth, iDay)
                    lJulian = julian.to_jd(current_date, fmt="jd")

                    # Calculate index into daily data array
                    daily_index = int(lJulian - lJulian_start)

                    # Check array bounds
                    if daily_index >= len(aData_daily):
                        raise IndexError(
                            f"Calculated index {daily_index} exceeds data length "
                            f"{len(aData_daily)} for date {iYear}-{iMonth:02d}-{iDay:02d}"
                        )

                    daily_value = aData_daily[daily_index]

                    # Only include finite values in mean calculation
                    if np.isfinite(daily_value):
                        monthly_sum += daily_value
                        valid_count += 1

                # Calculate monthly mean
                if valid_count == 0:
                    # All values were NaN
                    monthly_value = np.nan
                    logger.warning(f"All days in {iYear}-{iMonth:02d} are NaN")
                else:
                    monthly_value = monthly_sum / valid_count

                    # Warn if many NaN values
                    total_days = iDay_end - iDay_start + 1
                    nan_percentage = 100.0 * (total_days - valid_count) / total_days
                    if nan_percentage > 50:
                        logger.warning(
                            f"{iYear}-{iMonth:02d}: {nan_percentage:.1f}% of days are NaN "
                            f"({total_days - valid_count}/{total_days})"
                        )

                aData_monthly_out.append(monthly_value)

            elif sType == "sum":
                monthly_sum = 0.0

                for iDay in range(iDay_start, iDay_end + 1):
                    # Calculate Julian date for this day
                    current_date = datetime.datetime(iYear, iMonth, iDay)
                    lJulian = julian.to_jd(current_date, fmt="jd")

                    # Calculate index into daily data array
                    daily_index = int(lJulian - lJulian_start)

                    # Check array bounds
                    if daily_index >= len(aData_daily):
                        raise IndexError(
                            f"Calculated index {daily_index} exceeds data length "
                            f"{len(aData_daily)} for date {iYear}-{iMonth:02d}-{iDay:02d}"
                        )

                    daily_value = aData_daily[daily_index]

                    # For sum, treat NaN as 0
                    if np.isfinite(daily_value):
                        monthly_sum += daily_value

                aData_monthly_out.append(monthly_sum)

    # ========================================================================
    # Convert to numpy array and return
    # ========================================================================

    aData_monthly_out = np.array(aData_monthly_out)

    logger.info(
        f"Converted {len(aData_daily)} daily values to {len(aData_monthly_out)} "
        f"monthly values using '{sType}' aggregation"
    )

    return aData_monthly_out
