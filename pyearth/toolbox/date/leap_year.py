def leap_year(iYear_in):
    """
    Check whether this year is a leap year or not.

    Args:
        iYear_in (int): The year

    Returns:
        bool:  True or False
    """
    if iYear_in % 400 == 0:
        return True
    if iYear_in % 100 == 0:
        return False
    if iYear_in % 4 == 0:
        return True
    else:
        return False