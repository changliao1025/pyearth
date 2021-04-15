from jdcal import gcal2jd, jd2gcal

from pyearth.toolbox.date.leap_year import leap_year

def day_in_month(iYear_in, iMonth_in, iFlag_leap_year_in = None):
    lJulian_start = gcal2jd(iYear_in, iMonth_in, 1)

    if iMonth_in < 12:
        lJulian_end =  gcal2jd(iYear_in, iMonth_in+1, 1)
    else:
        lJulian_end =  gcal2jd(iYear_in+1, 1, 1)

    if iFlag_leap_year_in == 0:
        if leap_year(iYear_in):
            if iMonth_in == 2:
                dayinmon = 28
            else:
                dayinmon = int (lJulian_end[1] - lJulian_start[1] ) 
        else:

            dayinmon = int (lJulian_end[1] - lJulian_start[1] ) 
    else:


        dayinmon = int (lJulian_end[1] - lJulian_start[1] ) 
    return dayinmon