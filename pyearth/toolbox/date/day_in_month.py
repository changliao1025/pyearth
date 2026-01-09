"""
Days in Month Calculation

This module provides functionality for calculating the number of days in a given
month, with support for leap year handling. It uses Python's standard calendar
module for accurate date calculations across all calendar systems.

Main Functions
--------------
day_in_month : Calculate the number of days in a specific month

Key Features
------------
- Accurate day count for any month and year
- Automatic leap year detection
- Optional manual leap year override
- Input validation for month and year ranges
- Efficient calculation using standard library
- Comprehensive error handling

Use Cases
---------
1. **Date Range Calculations**: Determine month boundaries for data processing
2. **Calendar Applications**: Generate monthly views with correct day counts
3. **Time Series Analysis**: Calculate monthly intervals for scientific data
4. **Data Validation**: Verify date inputs fall within valid month ranges
5. **Hydrological Modeling**: Calculate monthly water balance with varying days
6. **Climate Analysis**: Process monthly climate data with correct day counts

Technical Details
-----------------
The module uses Python's calendar.monthrange() function for standard calculations,
which handles all leap year rules correctly:
- Years divisible by 4 are leap years
- Except years divisible by 100 are not leap years
- Except years divisible by 400 are leap years

For example:
- 2000: Leap year (divisible by 400)
- 1900: Not a leap year (divisible by 100 but not 400)
- 2024: Leap year (divisible by 4)
- 2023: Not a leap year

The iFlag_leap_year_in parameter allows manual override for specialized calendars
or hypothetical scenarios where standard leap year rules don't apply.

Performance Characteristics
---------------------------
- Time Complexity: O(1) - direct calculation
- Space Complexity: O(1) - minimal memory usage
- Very fast: uses optimized standard library functions

Dependencies
------------
- calendar: Standard library for calendar calculations
- logging: For error reporting and debugging

See Also
--------
- calendar.monthrange: Standard library function for day count
- calendar.isleap: Check if a year is a leap year
- datetime: Date and time manipulation
"""

import calendar
import logging
from typing import Optional

# Configure logging
logger = logging.getLogger(__name__)


def day_in_month(
    iYear_in: int, iMonth_in: int, iFlag_leap_year_in: Optional[int] = None
) -> int:
    """
    Calculate the number of days in a specific month.

    This function returns the number of days in a given month for a specified year,
    with automatic leap year detection. An optional parameter allows manual override
    of leap year status for specialized calendar systems or hypothetical scenarios.

    The function uses Python's calendar module for accurate calculations that properly
    handle all leap year rules including the Gregorian calendar exceptions.

    Parameters
    ----------
    iYear_in : int
        The year (e.g., 2024). Valid range: 1 to 9999.
        Used to determine leap year status for February.
    iMonth_in : int
        The month number (1-12), where:
        - 1 = January
        - 2 = February
        - 3 = March
        - ...
        - 12 = December
    iFlag_leap_year_in : Optional[int], optional
        Manual override for leap year status. If None (default), leap year is
        automatically detected using standard rules. If specified:
        - 0: Force non-leap year (February has 28 days)
        - 1: Force leap year (February has 29 days)
        Only affects February; other months are unaffected.

    Returns
    -------
    int
        The number of days in the specified month.
        - January: 31 days
        - February: 28 or 29 days (depends on leap year)
        - March: 31 days
        - April: 30 days
        - May: 31 days
        - June: 30 days
        - July: 31 days
        - August: 31 days
        - September: 30 days
        - October: 31 days
        - November: 30 days
        - December: 31 days

    Raises
    ------
    ValueError
        If iMonth_in is not in range 1-12.
        If iYear_in is not in valid range (1-9999).
        If iFlag_leap_year_in is not None, 0, or 1.

    Notes
    -----
    1. **Leap Year Rules**: The function follows Gregorian calendar leap year rules:
       - Years divisible by 4 are leap years
       - EXCEPT years divisible by 100 are not leap years
       - EXCEPT years divisible by 400 are leap years
       Examples: 2000 (leap), 1900 (not leap), 2024 (leap), 2100 (not leap)

    2. **Manual Override**: The iFlag_leap_year_in parameter allows overriding
       automatic leap year detection. This is useful for:
       - Hypothetical calendar scenarios
       - Non-Gregorian calendar systems
       - Testing and validation
       Only affects February; other months ignore this parameter.

    3. **Validation**: Input validation ensures month is 1-12 and year is in valid
       range. Invalid inputs raise ValueError with descriptive messages.

    4. **Month Day Counts**: Standard month lengths:
       - 31 days: Jan, Mar, May, Jul, Aug, Oct, Dec (7 months)
       - 30 days: Apr, Jun, Sep, Nov (4 months)
       - 28/29 days: Feb (1 month, varies by leap year)

    5. **Efficiency**: Uses calendar.monthrange() from standard library, which is
       highly optimized. No need for Julian date conversions or complex calculations.

    6. **Range Limitations**: Year must be in range 1-9999 due to datetime module
       limitations. This covers all practical use cases for modern applications.

    7. **Return Type**: Always returns an integer representing exact day count.
       Never returns None or raises exceptions for valid inputs.

    Examples
    --------
    Get days in February for a leap year (automatic detection):

    >>> days = day_in_month(2024, 2)
    >>> print(days)
    29

    Get days in February for a non-leap year:

    >>> days = day_in_month(2023, 2)
    >>> print(days)
    28

    Get days in a 31-day month:

    >>> days = day_in_month(2024, 1)  # January
    >>> print(days)
    31

    Get days in a 30-day month:

    >>> days = day_in_month(2024, 4)  # April
    >>> print(days)
    30

    Override leap year status (force non-leap year):

    >>> days = day_in_month(2024, 2, iFlag_leap_year_in=0)
    >>> print(days)
    28
    # Note: 2024 is actually a leap year, but override forces 28 days

    Override leap year status (force leap year):

    >>> days = day_in_month(2023, 2, iFlag_leap_year_in=1)
    >>> print(days)
    29
    # Note: 2023 is not a leap year, but override forces 29 days

    Verify edge case years:

    >>> # Year 2000: Leap year (divisible by 400)
    >>> print(day_in_month(2000, 2))
    29
    >>> # Year 1900: Not leap year (divisible by 100 but not 400)
    >>> print(day_in_month(1900, 2))
    28

    See Also
    --------
    calendar.monthrange : Get weekday and number of days for a month
    calendar.isleap : Check if a year is a leap year
    """
    # Validate month range
    if not isinstance(iMonth_in, int) or iMonth_in < 1 or iMonth_in > 12:
        logger.error(f"Invalid month: {iMonth_in}. Must be integer between 1 and 12.")
        raise ValueError(f"Month must be between 1 and 12, got {iMonth_in}")

    # Validate year range
    if not isinstance(iYear_in, int) or iYear_in < 1 or iYear_in > 9999:
        logger.error(f"Invalid year: {iYear_in}. Must be integer between 1 and 9999.")
        raise ValueError(f"Year must be between 1 and 9999, got {iYear_in}")

    # Validate leap year flag if provided
    if iFlag_leap_year_in is not None:
        if iFlag_leap_year_in not in [0, 1]:
            logger.error(
                f"Invalid leap year flag: {iFlag_leap_year_in}. Must be 0, 1, or None."
            )
            raise ValueError(
                f"Leap year flag must be 0, 1, or None, got {iFlag_leap_year_in}"
            )

    # Get number of days in month using standard library
    # calendar.monthrange returns (weekday, num_days)
    _, nDays = calendar.monthrange(iYear_in, iMonth_in)

    # Apply manual leap year override if specified and month is February
    if iFlag_leap_year_in is not None and iMonth_in == 2:
        if iFlag_leap_year_in == 0:
            nDays = 28
            logger.debug(
                f"Manual override: February {iYear_in} forced to 28 days (non-leap)"
            )
        else:  # iFlag_leap_year_in == 1
            nDays = 29
            logger.debug(
                f"Manual override: February {iYear_in} forced to 29 days (leap)"
            )

    logger.debug(f"Year {iYear_in}, Month {iMonth_in}: {nDays} days")

    return nDays
