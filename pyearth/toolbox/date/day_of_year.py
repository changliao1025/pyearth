import datetime
import pyearth.toolbox.date.julian as julian
def day_of_year(iYear_in, iMonth_in, iDay_in):
    """
    Calculate the day in a year

    Args:
        iYear_in (int): The year
        iMonth_in (int): The month
        iDay_in (int): The day

    Returns:
        int: The day in a year
    """

    if iYear_in is None:
        print("Parameter 'iYear_in' is undefined ")
        return -1
    if iMonth_in is None:
        print("Parameter 'iMonth_in' is undefined ")
        return -1
    else:
       if iMonth_in < 1 or iMonth_in > 12:
          print("Parameter 'iMonth_in' is not within correct range ")
          return -1
    if iDay_in is None:
        print("Parameter 'iDay_in' is undefined ")
        return -1
    else:
       if iDay_in < 1 or iDay_in > 31:
          print('Parameter iDay_in is not within correct range')
          return -1

    dummy1 = datetime.datetime(iYear_in, 1, 1)
    lJulian_start = julian.to_jd(dummy1, fmt='jd')

    dummy1 = datetime.datetime(iYear_in, iMonth_in, iDay_in)
    lJulian_end = julian.to_jd(dummy1, fmt='jd')

    dayofyear = int (lJulian_end[1] - lJulian_start[1] ) + 1

    return dayofyear
