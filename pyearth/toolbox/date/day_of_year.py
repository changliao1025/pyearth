#import julian
#import datetime

from jdcal import gcal2jd, jd2gcal

def day_of_year(iYear_in, iMonth_in, iDay_in):

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

          
    #first_day_of_year = datetime.datetime(iYear_in, 1, 1)
    #current_day_of_year = datetime.datetime(iYear_in, iMonth_in, iDay_in)

    #lJulian_start = julian.to_jd(first_day_of_year, fmt='jd')
    #lJulian_end = julian.to_jd(current_day_of_year, fmt='jd')

    lJulian_start = gcal2jd(iYear_in, 1, 1)
    lJulian_end =  gcal2jd(iYear_in, iMonth_in, iDay_in)


    dayofyear = int (lJulian_end[1] - lJulian_start[1] ) + 1
    return dayofyear
