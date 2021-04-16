from jdcal import gcal2jd, jd2gcal


def day_in_month(iYear_in, iMonth_in):
    lJulian_start = gcal2jd(iYear_in, iMonth_in, 1)

    """Calculate how many days in a month"""

    if iMonth_in < 12:
        lJulian_end =  gcal2jd(iYear_in, iMonth_in+1, 1)
    else:
        lJulian_end =  gcal2jd(iYear_in+1, 1, 1)

    
    
   
    dayinmon = int (lJulian_end[1] - lJulian_start[1] ) 
  

    return dayinmon