

def retrieve_usgs_site_information_nwis(sSiteId):
    try:
        import requests
    except ImportError as e:
        raise ImportError(
            "The package 'requests' is required for this function to run.") from e
    
    squaremile2squarekm = 2.58999
    
    nwis_site_url = f"https://waterservices.usgs.gov/nwis/site/?format=rdb&siteOutput=Expanded&sites={sSiteId}"
    response = requests.get(nwis_site_url)
    pData = response.text       
    dLon = None
    dLat = None
    dDrainage = None
    
    for line in pData.splitlines():
        if line.startswith("#"):
            #print(line)
            continue
        data = line.split("\t")
        if len(data) >= 4:
            #print(data)
            if data[0] == 'agency_cd':
                index_lon = data.index('dec_long_va')
                index_lat = data.index('dec_lat_va')
                index_drai = data.index('drain_area_va')      
            
            if data[0] == 'USGS' and data[1] == sSiteId:
                dLon = float(data[index_lon])
                dLat = float(data[index_lat])
                dDrainage= float(data[index_drai]) * squaremile2squarekm
    
    return dLon, dLat, dDrainage