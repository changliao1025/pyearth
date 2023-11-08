import numpy as np
from datetime import datetime

from pyearth.toolbox.date.day_in_month import day_in_month

def retrieve_nwis_discharge(sSite, 
                            iYear_start, 
                            iMonth_start, 
                            iDay_start,  
                            iYear_end, 
                            iMonth_end, 
                            iDay_end):
    try:
        import requests
    except ImportError as e:
        raise ImportError(
            "The package 'requests' is required for this function to run.") from e
    
    nwis_dv_url = 'https://waterservices.usgs.gov/nwis/dv/'

    aDate= list()
    sYear_start = "{:04d}".format(iYear_start)
    sMonth_start = "{:02d}".format(iMonth_start)
    sDay_start = "{:02d}".format(iDay_start)
    sYear_end = "{:04d}".format(iYear_end)
    sMonth_end = "{:02d}".format(iMonth_end)
    sDay_end = "{:02d}".format(iDay_end)

    for iYear in range(iYear_start,iYear_end+1):
        for iMonth in range(1,13):
            iDay_end = day_in_month(iYear,iMonth)
            for iDay in range(1,iDay_end+1):            
                aDate.append(datetime(iYear,iMonth,iDay))
    nDay = len(aDate)
    aDate = np.array(aDate)

    start_date = sYear_start + '-' + sMonth_start + '-' + sDay_start
    end_date = sYear_end + '-' + sMonth_end + '-' + sDay_end

    params = {
        'format': 'json',
        'sites': sSite,
        'startDT': start_date,
        'endDT': end_date,
        'parameterCd': '00060' # USGS code for streamflow

    }
    # Send API request
    aDischarge = np.full(nDay, np.nan, dtype= float)
    response = requests.get(nwis_dv_url, params=params)
    unit = None
    # Check if request was successful
    if response.status_code == 200:
        # Parse JSON response
        data = response.json()    
        # Extract streamflow values
        time_series = data['value']['timeSeries'][0]['values'][0]['value']
        unit = data['value']['timeSeries'][0]['variable']['unit']['unitCode']
        print("Unit of measurement:", unit)    
        # Iterate over streamflow values
        for value in time_series:
            date = value['dateTime']
            #take the substring
            date_string = date[0:10]
            # Convert the string to datetime
            iYear = int(date_string[0:4])
            iMonth = int(date_string[5:7])
            iDay = int(date_string[8:10])
            datetime_object = datetime(iYear,iMonth,iDay)
            index = np.where( aDate == datetime_object)
            aDischarge[index] = value['value']  
    else:
        print('Error:', response.status_code)

    return aDischarge, unit