import numpy as np
from jdcal import gcal2jd, jd2gcal
from pyearth.toolbox.date.day_in_month import day_in_month
def convert_time_series_daily_to_monthly(aData_daily_in,\
    iYear_start_in, iMonth_start_in, iDay_start_in, \
      iYear_end_in, iMonth_end_in, iDay_end_in , sType_in = None  ):
    
    if sType_in is not None:
        sType = sType_in
    else:
        sType = 'mean'

    #check the length of data

    aData_daily = np.arra(aData_daily_in)
    aData_monthly_out = list()

    lJulian_start = gcal2jd(iYear_start_in, iMonth_start_in, iDay_start_in)
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
                    lJulian=gcal2jd(iYear, iMonth, iDay)
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
                        lJulian=gcal2jd(iYear, iMonth, iDay)
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
