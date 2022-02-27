import datetime
import pyearth.toolbox.date.julian as julian
def day_in_month(iYear_in, iMonth_in):
    """
    Calculate how many days in a month

    Args:
        iYear_in (int): The year
        iMonth_in (int): The month

    Returns:
        int: The total day in the month
    """
    
    dummy1 = datetime.datetime(iYear_in, iMonth_in, 1)
    lJulian_start = julian.to_jd(dummy1, fmt='jd')    

    if iMonth_in < 12:
        dummy1 = datetime.datetime(iYear_in, iMonth_in+1, 1)
        lJulian_end = julian.to_jd(dummy1, fmt='jd')
        
    else:
        dummy1 = datetime.datetime(iYear_in+1, iMonth_in, 1)
        lJulian_end = julian.to_jd(dummy1, fmt='jd')
   
    dayinmon = int (lJulian_end - lJulian_start ) 
  

    return dayinmon