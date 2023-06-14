import requests

def retrieve_usgs_site_information_nldi(sSiteID):
    
    nldi_site_url = f"https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/USGS-{sSiteID}"
    response = requests.get(nldi_site_url)   
    pData = response.json()    
    aLocation  = pData['features'][0]['geometry']['coordinates']
    dLon = aLocation[0]
    dLat = aLocation[1]
    
    return dLon, dLat