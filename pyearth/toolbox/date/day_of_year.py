"""
Day of Year Calculation

This module provides functionality for calculating the day of year (ordinal day)
for a given date. The day of year is a number between 1 and 366 representing the
sequential day count from January 1st of the year.

Main Functions
--------------
day_of_year : Calculate the ordinal day number for a specific date

Key Features
------------
- Accurate day of year calculation for any valid date
- Automatic leap year handling (February 29 = day 60)
- Comprehensive input validation with informative error messages
- Efficient calculation using Python's standard datetime module
- Type hints for better code clarity
- Professional logging for debugging

Use Cases
---------
1. **Time Series Analysis**: Convert dates to sequential day numbers for analysis
2. **Climate Data Processing**: Organize meteorological data by day of year
3. **Seasonal Analysis**: Group data by day of year across multiple years
4. **Data Indexing**: Use ordinal day as an index for daily datasets
5. **Phenology Studies**: Track seasonal events using day of year
6. **Agricultural Planning**: Schedule activities based on day of year
7. **Solar Position**: Calculate sun position requires day of year

Technical Details
-----------------
The day of year (also called ordinal day or Julian day) is the sequential count
of days from January 1st of a given year:
- January 1 = Day 1
- December 31 = Day 365 (non-leap years) or Day 366 (leap years)
- February 29 = Day 60 (leap years only)

This module uses Python's datetime.timetuple().tm_yday for calculation, which
automatically handles all calendar complexities including leap years.

Leap Year Rules (Gregorian Calendar):
- Years divisible by 4 are leap years (366 days)
- EXCEPT years divisible by 100 are not leap years
- EXCEPT years divisible by 400 are leap years

Examples:
- 2024: Leap year (366 days, Feb 29 exists)
- 2023: Regular year (365 days)
- 2000: Leap year (divisible by 400)
- 1900: Not a leap year (divisible by 100 but not 400)

Performance Characteristics
---------------------------
- Time Complexity: O(1) - direct calculation
- Space Complexity: O(1) - minimal memory usage
- Very fast: uses optimized datetime module

Dependencies
------------
- datetime: Standard library for date/time handling
- logging: For error reporting and debugging

See Also
--------
- datetime.date.timetuple: Get time tuple including day of year
- datetime.date.toordinal: Get proleptic Gregorian ordinal
- day_in_month: Calculate number of days in a specific month
"""

import datetime
import logging
from typing import Union

# Configure logging
logger = logging.getLogger(__name__)


def day_of_year(iYear_in: int, iMonth_in: int, iDay_in: int) -> int:
    """
    Calculate the day of year (ordinal day) for a specific date.

    This function returns the sequential day number from January 1st for a given
    date. The result ranges from 1 (January 1) to 365 (December 31 in non-leap
    years) or 366 (December 31 in leap years).

    The calculation automatically handles leap years according to Gregorian
    calendar rules, ensuring February 29 is correctly handled as day 60 in
    leap years.

    Parameters
    ----------
    iYear_in : int
        The year (e.g., 2024). Valid range: 1 to 9999.
        Required to determine leap year status.
    iMonth_in : int
        The month number (1-12), where:
        - 1 = January
        - 2 = February
        - ...
        - 12 = December
    iDay_in : int
        The day of the month (1-31). Must be valid for the given month.
        For example, February 30 is invalid and will raise ValueError.

    Returns
    -------
    int
        The day of year (ordinal day) ranging from:
        - 1: January 1st
        - 32: February 1st (non-leap year) or February 1st (leap year)
        - 60: February 29th (leap years only)
        - 365: December 31st (non-leap years)
        - 366: December 31st (leap years)

    Raises
    ------
    ValueError
        If any parameter is None, out of valid range, or represents an invalid date.
        Examples of invalid dates:
        - February 30 (month only has 28/29 days)
        - April 31 (month only has 30 days)
        - February 29, 2023 (2023 is not a leap year)
    TypeError
        If parameters are not integers.

    Notes
    -----
    1. **Day of Year Range**: Results range from 1 to 365 (non-leap years) or
       1 to 366 (leap years). There is no day 0.

    2. **Leap Year Handling**: The function automatically handles leap years using
       Python's datetime module, which follows Gregorian calendar rules:
       - Years divisible by 4: Leap years (366 days)
       - Except divisible by 100: Not leap years
       - Except divisible by 400: Leap years

    3. **Date Validation**: The datetime module validates that the date is valid
       for the given year. For example, February 29, 2023 raises ValueError because
       2023 is not a leap year.

    4. **Common Day Numbers**:
       - January 1: Day 1
       - March 1: Day 60 (leap) or Day 59 (non-leap)
       - July 1: Day 182 (leap) or Day 181 (non-leap)
       - December 31: Day 366 (leap) or Day 365 (non-leap)

    5. **Efficient Implementation**: Uses datetime.timetuple().tm_yday which is
       highly optimized in the CPython implementation. Much faster than manual
       calculation or Julian date conversions.

    6. **Year Range**: Limited to 1-9999 by Python's datetime module. This covers
       all practical modern applications.

    7. **Comparison with Ordinal**: This function returns day within year (1-366),
       not proleptic Gregorian ordinal. For ordinal, use datetime.date.toordinal().

    Examples
    --------
    Calculate day of year for New Year's Day:

    >>> doy = day_of_year(2024, 1, 1)
    >>> print(doy)
    1

    Calculate day of year for February 29 (leap year):

    >>> doy = day_of_year(2024, 2, 29)
    >>> print(doy)
    60

    Calculate day of year for March 1 in leap year:

    >>> doy = day_of_year(2024, 3, 1)
    >>> print(doy)
    61
    # Note: Day 61 because Feb 29 exists in 2024

    Calculate day of year for March 1 in non-leap year:

    >>> doy = day_of_year(2023, 3, 1)
    >>> print(doy)
    60
    # Note: Day 60 because Feb 29 doesn't exist in 2023

    Calculate day of year for December 31 (leap year):

    >>> doy = day_of_year(2024, 12, 31)
    >>> print(doy)
    366

    Calculate day of year for December 31 (non-leap year):

    >>> doy = day_of_year(2023, 12, 31)
    >>> print(doy)
    365

    Calculate day of year for Independence Day:

    >>> doy = day_of_year(2024, 7, 4)
    >>> print(doy)
    186

    Verify leap year edge cases:

    >>> # Year 2000: Leap year (divisible by 400)
    >>> print(day_of_year(2000, 2, 29))
    60
    >>> # Year 1900: Not leap year (divisible by 100 but not 400)
    >>> # day_of_year(1900, 2, 29)  # Would raise ValueError
    >>> print(day_of_year(1900, 3, 1))
    60

    See Also
    --------
    day_in_month : Calculate number of days in a specific month
    datetime.date.timetuple : Get time tuple with day of year
    datetime.date.toordinal : Get proleptic Gregorian ordinal
    """
    # Validate inputs are not None
    if iYear_in is None:
        logger.error("Parameter 'iYear_in' is None")
        raise ValueError("Year cannot be None")

    if iMonth_in is None:
        logger.error("Parameter 'iMonth_in' is None")
        raise ValueError("Month cannot be None")

    if iDay_in is None:
        logger.error("Parameter 'iDay_in' is None")
        raise ValueError("Day cannot be None")

    # Validate inputs are integers
    if not isinstance(iYear_in, int):
        logger.error(f"Year must be integer, got {type(iYear_in).__name__}")
        raise TypeError(f"Year must be integer, got {type(iYear_in).__name__}")

    if not isinstance(iMonth_in, int):
        logger.error(f"Month must be integer, got {type(iMonth_in).__name__}")
        raise TypeError(f"Month must be integer, got {type(iMonth_in).__name__}")

    if not isinstance(iDay_in, int):
        logger.error(f"Day must be integer, got {type(iDay_in).__name__}")
        raise TypeError(f"Day must be integer, got {type(iDay_in).__name__}")

    # Validate year range
    if iYear_in < 1 or iYear_in > 9999:
        logger.error(f"Invalid year: {iYear_in}. Must be between 1 and 9999.")
        raise ValueError(f"Year must be between 1 and 9999, got {iYear_in}")

    # Validate month range
    if iMonth_in < 1 or iMonth_in > 12:
        logger.error(f"Invalid month: {iMonth_in}. Must be between 1 and 12.")
        raise ValueError(f"Month must be between 1 and 12, got {iMonth_in}")

    # Validate day range (basic check, datetime will do detailed validation)
    if iDay_in < 1 or iDay_in > 31:
        logger.error(f"Invalid day: {iDay_in}. Must be between 1 and 31.")
        raise ValueError(f"Day must be between 1 and 31, got {iDay_in}")

    # Create date object and calculate day of year
    # This will raise ValueError if date is invalid (e.g., Feb 30, Feb 29 in non-leap year)
    try:
        date_obj = datetime.date(iYear_in, iMonth_in, iDay_in)
        nDay_of_year = date_obj.timetuple().tm_yday

        logger.debug(
            f"Date {iYear_in}-{iMonth_in:02d}-{iDay_in:02d} is day {nDay_of_year} of year"
        )

        return nDay_of_year

    except ValueError as e:
        logger.error(f"Invalid date {iYear_in}-{iMonth_in}-{iDay_in}: {e}")
        raise ValueError(
            f"Invalid date: {iYear_in}-{iMonth_in:02d}-{iDay_in:02d}. {str(e)}"
        )
