#write the header for ENVI raster file from the dictionary data type
#https://www.harrisgeospatial.com/docs/ENVIGridDefinition.html
def envi_write_header(sFilename_in, aHeader):
    pFile = open(sFilename_in, 'w')
    sLine = 'ENVI' + '\n'
    pFile.write(sLine)
    sLine = 'description = ' + aHeader['sFilename'] + '\n'
    pFile.write(sLine)
    sLine = 'samples = '  + aHeader['ncolumn'] + '\n'
    pFile.write(sLine)
    sLine = 'lines = ' + aHeader['nrow'] + '\n'
    pFile.write(sLine)
    sLine = 'bands = ' + aHeader['nband'] + '\n'
    pFile.write(sLine)
    sLine = 'header offset = ' + aHeader['offset'] + '\n'
    pFile.write(sLine)
    sLine = 'data type = ' + aHeader['data_type'] + '\n'
    pFile.write(sLine)
    sLine = 'interleave = ' + aHeader['bsq'] + '\n'
    pFile.write(sLine)
    sLine = 'sensor type = Unknown' + '\n'
    pFile.write(sLine)
    sLine = 'byte order = ' + aHeader['byte_order'] + '\n'
    pFile.write(sLine)
    sLine = 'data ignore value = ' + aHeader['missing_value'] + '\n'
    pFile.write(sLine)
    sLine = 'map info = {Geographic Lat/Lon, 1.000, 1.000, ' \
        + aHeader['ULlon'] + ', ' + aHeader['ULlat'] + ', ' \
        + aHeader['pixelSize'] + ', ' + aHeader['pixelSize'] + ', ' \
        + 'WGS-84, units = Degrees}' + '\n'
    pFile.write(sLine)
    sLine = 'wavelength units = Unknown ' + '\n'
    pFile.write(sLine)
    pFile.close()