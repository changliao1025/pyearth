import numpy as np
import datetime
import pyearth.toolbox.date.julian as julian
from pyearth.toolbox.date.day_in_month import day_in_month

def convert_time_series_daily_to_monthly(aData_daily_in, iYear_start_in, iMonth_start_in, iDay_start_in, iYear_end_in, iMonth_end_in, iDay_end_in , sType_in = None  ):
    """
    Convert a daily time series data into a monthly data

    Args:
        aData_daily_in (numpy.array): The daily time series data
        iYear_start_in (int): The start year
        iMonth_start_in (int): The start month
        iDay_start_in (int): The start day
        iYear_end_in (int): The end year
        iMonth_end_in (int): The end month
        iDay_end_in (int): The end day
        sType_in (string, optional): Method for aggregation. Defaults to None.

    Returns:
        numpy.array: aData_monthly_out
    """
    if sType_in is not None:
        sType = sType_in
    else:
        sType = 'mean'

    #check the length of data

    aData_daily = np.array(aData_daily_in)
    aData_monthly_out = list()

    dummy1 = datetime.datetime(iYear_start_in, iMonth_start_in, iDay_start_in)
    #dummy2 = datetime.datetime(iYear_end, iMonth_end, iDay_end)
    lJulian_start = julian.to_jd(dummy1, fmt='jd')
    #julian2 = julian.to_jd(dummy2, fmt='jd')
    #the old package
    #lJulian_start = gcal2jd(iYear_start_in, iMonth_start_in, iDay_start_in)
    for iYear in range( iYear_start_in, iYear_end_in +1):
        if iYear == iYear_start_in:
            iMonth_start = iMonth_start_in
        else:
            iMonth_start = 1

        if iYear == iYear_end_in :
            iMonth_end = iMonth_end_in
        else:
            iMonth_end = 12

        for iMonth in range(iMonth_start, iMonth_end+1):
            #get all the data from this year-month

            if sType == "mean":

                if iYear == iYear_start_in and iMonth ==iMonth_start_in:
                    iDay_start = iDay_start_in
                else: 
                    iDay_start = 1
                
                if iYear == iYear_end_in and iMonth ==iMonth_end_in:
                    iDay_end = iDay_end_in
                else: 
                    iDay_end = day_in_month(iYear, iMonth)
                
                dummy = 0.0
                count = 0
                for iDay in range(iDay_start, iDay_end+1):
                    dummy2 = datetime.datetime(iYear, iMonth, iDay)
                    lJulian = julian.to_jd(dummy2, fmt='jd')
                    #lJulian=gcal2jd(iYear, iMonth, iDay)
                    dummy_index = int(lJulian-lJulian_start)
                    dummy = dummy + aData_daily[dummy_index]
                    count = count + 1
                    pass

                aData_monthly_out.append( dummy/count )

                pass
            else:
                if sType == 'sum':
                    if iYear == iYear_start_in and iMonth ==iMonth_start_in:
                        iDay_start = iDay_start_in
                    else: 
                        iDay_start = 1

                    if iYear == iYear_end_in and iMonth ==iMonth_end_in:
                        iDay_end = iDay_end_in
                    else: 
                        iDay_end = day_in_month(iYear, iMonth)

                    dummy = 0.0
                    
                    for iDay in range(iDay_start, iDay_end+1):
                        #lJulian=gcal2jd(iYear, iMonth, iDay)
                        dummy2 = datetime.datetime(iYear, iMonth, iDay)
                        lJulian = julian.to_jd(dummy2, fmt='jd')
                        dummy_index = int(lJulian-lJulian_start)
                        dummy = dummy + aData_daily[dummy_index]
                        pass

                    aData_monthly_out.append( dummy )
                    pass
                else:
                    pass




            pass

    aData_monthly_out = np.array(aData_monthly_out)
    return aData_monthly_out
