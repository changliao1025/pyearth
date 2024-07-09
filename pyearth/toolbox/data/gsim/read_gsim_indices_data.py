import numpy as np
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from pyearth.toolbox.date.day_in_month import day_in_month

def read_gsim_indices_data(sFilename_in, iYear_start_in=2000, iYear_end_in=2019):

    dummy = text_reader_string(sFilename_in, iSkipline_in=22, cDelimiter_in=',')

    aDate = dummy[:, 0]
    aDischarge_mean = np.array([float(x.strip()) if x.strip() != 'NA' else np.nan for x in dummy[:, 1]])

    #build the monthly index

    nYear = iYear_end_in - iYear_start_in + 1

    nstress = 12 * nYear

    aDate_out = np.full(nstress, object)

    aDischarge_out = np.full(nstress, np.nan)

    i = 0

    for iYear in range(iYear_start_in, iYear_end_in + 1):
        for iMonth in range(1, 13):
            nday = day_in_month(iYear, iMonth)
            sDate = "{}-{}-{}".format(iYear, str(iMonth).zfill(2), str(nday).zfill(2))
            iIndex = np.where(aDate == sDate)
            if len(iIndex[0]) == 0:
                aDischarge_out[i] = np.nan
            else:
                aDischarge_out[i] = aDischarge_mean[iIndex]

            aDate_out[i] = sDate
            i = i + 1


    return aDate_out, aDischarge_out
