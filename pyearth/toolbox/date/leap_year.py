def leap_year(y):
    """
    Check whether this year is a leap year or not.
    """
    if y % 400 == 0:
        return True
    if y % 100 == 0:
        return False
    if y % 4 == 0:
        return True
    else:
        return False