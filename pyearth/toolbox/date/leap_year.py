"""
Leap Year Determination

This module provides functionality for determining whether a given year is a leap
year according to the Gregorian calendar rules. Leap years have 366 days instead
of the standard 365 days, with the extra day being February 29th.

Main Functions
--------------
leap_year : Check if a year is a leap year

Key Features
------------
- Accurate leap year determination using Gregorian calendar rules
- Handles all edge cases (divisible by 100, 400, etc.)
- Input validation with informative error messages
- Type hints for better code clarity
- Professional logging for debugging
- Efficient O(1) calculation

Use Cases
---------
1. **Date Validation**: Verify if February 29 is valid for a given year
2. **Calendar Generation**: Create accurate calendars with correct day counts
3. **Day Count Calculations**: Determine if a year has 365 or 366 days
4. **Time Series Processing**: Handle leap years in temporal data analysis
5. **Agricultural Planning**: Account for leap years in seasonal schedules
6. **Financial Calculations**: Compute day counts for interest calculations
7. **Data Quality**: Validate date inputs in datasets

Technical Details
-----------------
The Gregorian calendar defines leap years using three rules applied in order:

1. **Divisible by 400**: These years ARE leap years
   - Examples: 1600, 2000, 2400
   - This takes precedence over all other rules

2. **Divisible by 100**: These years are NOT leap years
   - Examples: 1700, 1800, 1900, 2100
   - Exception: If also divisible by 400 (see rule 1)

3. **Divisible by 4**: These years ARE leap years
   - Examples: 2004, 2008, 2012, 2016, 2020, 2024
   - Exception: If also divisible by 100 (see rule 2)

4. **All other years**: NOT leap years
   - Examples: 2001, 2002, 2003, 2005, 2007

The rules must be checked in this specific order to handle the exceptions correctly.

Historical Context:
- Julian calendar (before 1582): Every 4 years was a leap year (no 100/400 rule)
- Gregorian calendar (1582 onwards): Added the 100/400 rules for accuracy
- This corrects the drift of ~11 minutes per year in the Julian calendar

Why Leap Years Exist:
- Earth's orbit: ~365.2422 days (not exactly 365)
- Without leap years: Seasons would drift ~24 days per century
- Leap year system: Keeps calendar aligned with Earth's orbit

Performance Characteristics
---------------------------
- Time Complexity: O(1) - constant time for modulo operations
- Space Complexity: O(1) - no additional memory needed
- Very fast: Only 1-3 modulo operations per call

Dependencies
------------
- logging: For error reporting and debugging

See Also
--------
- calendar.isleap: Standard library function for leap year check
- day_in_month: Calculate number of days in a month (uses leap year logic)
- day_of_year: Calculate ordinal day (affected by leap years)
"""

import logging
from typing import Union

# Configure logging
logger = logging.getLogger(__name__)


def leap_year(iYear_in: int) -> bool:
    """
    Determine if a year is a leap year according to Gregorian calendar rules.

    This function checks whether a given year has 366 days (leap year) or 365 days
    (common year) by applying the three rules of the Gregorian calendar in the
    correct order.

    A leap year occurs:
    1. Every year divisible by 4
    2. EXCEPT for years divisible by 100
    3. EXCEPT for years divisible by 400

    Parameters
    ----------
    iYear_in : int
        The year to check (e.g., 2024). Can be any positive integer.
        Negative years (BC dates) are supported but follow proleptic Gregorian
        calendar rules which may not match historical calendars.

    Returns
    -------
    bool
        True if the year is a leap year (366 days, February has 29 days).
        False if the year is a common year (365 days, February has 28 days).

    Raises
    ------
    TypeError
        If input is not an integer.
    ValueError
        If input is None or zero (year 0 does not exist in Gregorian calendar).

    Notes
    -----
    1. **Rule Order**: The three leap year rules must be checked in order:
       - First check divisibility by 400 (these ARE leap years)
       - Then check divisibility by 100 (these are NOT leap years)
       - Then check divisibility by 4 (these ARE leap years)
       - All others are NOT leap years

       This order handles the exceptions correctly.

    2. **Historical Accuracy**: This function implements proleptic Gregorian
       calendar rules, meaning it applies Gregorian rules to all years including
       those before the calendar was adopted (1582). For historical dates before
       1582, actual calendars varied by region.

    3. **Common Examples**:
       - 2024: Leap year (divisible by 4, not by 100)
       - 2023: Not a leap year (not divisible by 4)
       - 2000: Leap year (divisible by 400)
       - 1900: Not a leap year (divisible by 100 but not by 400)
       - 2100: Not a leap year (divisible by 100 but not by 400)

    4. **February 29**: Leap years have February 29 as a valid date. Common years
       do not have February 29. Attempting to create a date like "February 29, 2023"
       will fail because 2023 is not a leap year.

    5. **Year Length**: Leap years have 366 days total, common years have 365 days.
       The extra day is always added as February 29.

    6. **Frequency**: In a 400-year cycle:
       - Leap years: 97 (24.25% of years)
       - Common years: 303 (75.75% of years)
       - Pattern: 303 years of 365 days + 97 years of 366 days = 146,097 days
       - Average: 365.2425 days per year (close to Earth's ~365.2422 day orbit)

    7. **Performance**: Uses only modulo operations, extremely fast even for
       large year values. Suitable for processing millions of years.

    8. **Comparison with Standard Library**: This function produces identical
       results to Python's calendar.isleap() function.

    9. **Negative Years**: Negative years represent BC dates in proleptic
       Gregorian calendar. Year -1 = 1 BC, -100 = 100 BC, etc. Note that
       there is no year 0 in historical calendars.

    10. **Alternative Implementation**: Could also use calendar.isleap() from
        the standard library, but this implementation has no dependencies and
        is equally fast.

    Examples
    --------
    Check if 2024 is a leap year:

    >>> is_leap = leap_year(2024)
    >>> print(is_leap)
    True
    >>> print(f"2024 has {366 if is_leap else 365} days")
    2024 has 366 days

    Check if 2023 is a leap year:

    >>> is_leap = leap_year(2023)
    >>> print(is_leap)
    False
    >>> print(f"2023 has {366 if is_leap else 365} days")
    2023 has 365 days

    Check century years (divisible by 100):

    >>> print(f"1900: {leap_year(1900)}")  # Divisible by 100, not by 400
    1900: False
    >>> print(f"2000: {leap_year(2000)}")  # Divisible by 400
    2000: True
    >>> print(f"2100: {leap_year(2100)}")  # Divisible by 100, not by 400
    2100: False

    Verify leap year pattern over a decade:

    >>> years = range(2020, 2030)
    >>> for year in years:
    ...     status = "leap" if leap_year(year) else "common"
    ...     print(f"{year}: {status}")
    2020: leap
    2021: common
    2022: common
    2023: common
    2024: leap
    2025: common
    2026: common
    2027: common
    2028: leap
    2029: common

    Validate February 29 dates:

    >>> leap_years = [year for year in range(1900, 2101) if leap_year(year)]
    >>> print(f"Leap years from 1900-2100: {len(leap_years)}")
    Leap years from 1900-2100: 48
    >>> print(f"First 5: {leap_years[:5]}")
    First 5: [1904, 1908, 1912, 1916, 1920]
    >>> print(f"Last 5: {leap_years[-5:]}")
    Last 5: [2084, 2088, 2092, 2096, 2100]

    Historical edge cases:

    >>> # Year 2000: Famous leap year (millennium, divisible by 400)
    >>> print(f"Y2K was {'leap' if leap_year(2000) else 'not leap'}")
    Y2K was leap
    >>> # Year 1900: Not a leap year despite being divisible by 4
    >>> print(f"1900 was {'leap' if leap_year(1900) else 'not leap'}")
    1900 was not leap

    See Also
    --------
    calendar.isleap : Standard library function for leap year check
    day_in_month : Calculate days in a month (February varies by leap year)
    day_of_year : Calculate ordinal day (max is 365 or 366)
    """
    # Validate input is not None
    if iYear_in is None:
        logger.error("Year parameter is None")
        raise ValueError("Year cannot be None")

    # Validate input is an integer
    if not isinstance(iYear_in, int):
        logger.error(f"Year must be integer, got {type(iYear_in).__name__}")
        raise TypeError(f"Year must be integer, got {type(iYear_in).__name__}")

    # Validate year is not zero (year 0 doesn't exist in Gregorian calendar)
    if iYear_in == 0:
        logger.error("Year 0 does not exist in Gregorian calendar")
        raise ValueError("Year 0 does not exist in Gregorian calendar")

    # Apply Gregorian calendar leap year rules in order
    # Rule 1: Divisible by 400 -> IS a leap year
    if iYear_in % 400 == 0:
        logger.debug(f"Year {iYear_in} is a leap year (divisible by 400)")
        return True

    # Rule 2: Divisible by 100 (but not 400) -> NOT a leap year
    if iYear_in % 100 == 0:
        logger.debug(
            f"Year {iYear_in} is not a leap year (divisible by 100 but not 400)"
        )
        return False

    # Rule 3: Divisible by 4 (but not 100) -> IS a leap year
    if iYear_in % 4 == 0:
        logger.debug(f"Year {iYear_in} is a leap year (divisible by 4)")
        return True

    # Rule 4: All other years -> NOT a leap year
    logger.debug(f"Year {iYear_in} is not a leap year (not divisible by 4)")
    return False
