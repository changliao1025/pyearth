import os
def envi_write_header(sFilename_in, aHeader_in):
    """
    write the header for ENVI raster file from the dictionary data type, 
    currently it only supports WGS84 and a dict is required to store all the information. 
    https://www.harrisgeospatial.com/docs/ENVIGridDefinition.html

    Parameters
    ----------
        sFilename_in (string): The filename
        aHeader (dict): A dict that contains all required information
    Returns
    -------

    See Also
    --------
    
    Examples
    --------
    """
    
    
    if os.path.exists(sFilename_in):
        pass
    else:
        print('The envi file does not exist!')
        return

    pFile = open(sFilename_in, 'w')
    sLine = 'ENVI' + '\n'
    pFile.write(sLine)
    sLine = 'description = ' + aHeader_in['sFilename'] + '\n'
    pFile.write(sLine)
    sLine = 'samples = '  + aHeader_in['ncolumn'] + '\n'
    pFile.write(sLine)
    sLine = 'lines = ' + aHeader_in['nrow'] + '\n'
    pFile.write(sLine)
    sLine = 'bands = ' + aHeader_in['nband'] + '\n'
    pFile.write(sLine)
    sLine = 'header offset = ' + aHeader_in['offset'] + '\n'
    pFile.write(sLine)
    sLine = 'data type = ' + aHeader_in['data_type'] + '\n'
    pFile.write(sLine)
    sLine = 'interleave = ' + aHeader_in['bsq'] + '\n'
    pFile.write(sLine)
    sLine = 'sensor type = Unknown' + '\n'
    pFile.write(sLine)
    sLine = 'byte order = ' + aHeader_in['byte_order'] + '\n'
    pFile.write(sLine)
    sLine = 'data ignore value = ' + aHeader_in['missing_value'] + '\n'
    pFile.write(sLine)

    #only one spatial reference is supported here.
    sLine = 'map info = {Geographic Lat/Lon, 1.000, 1.000, ' \
        + aHeader_in['ULlon'] + ', ' + aHeader_in['ULlat'] + ', ' \
        + aHeader_in['pixelSize'] + ', ' + aHeader_in['pixelSize'] + ', ' \
        + 'WGS-84, units = Degrees}' + '\n'
    pFile.write(sLine)

    sLine = 'wavelength units = Unknown ' + '\n'
    pFile.write(sLine)
    pFile.close()
