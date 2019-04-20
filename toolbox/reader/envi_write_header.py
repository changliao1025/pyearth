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
    sLine = 'bands = 1 ' + '\n'
    pFile.write(sLine)
    sLine = 'header offset = 0' + '\n'
    pFile.write(sLine)
    sLine = 'data type = 4' + '\n'
    pFile.write(sLine)
    sLine = 'interleave = bsq' + '\n'
    pFile.write(sLine)
    sLine = 'sensor type = Unknown' + '\n'
    pFile.write(sLine)
    sLine = 'byte order = 0' + '\n'
    pFile.write(sLine)
    sLine = 'data ignore value = -9999' + '\n'
    pFile.write(sLine)
    sLine = '{Geographic Lat/Lon, 1.000, 1.000, ' + aHeader['ULlon'] + ', ' + aHeader['ULlat'] + ', ' + aHeader['pixelSize'] + ', '+ aHeader['pixelSize'] + ', WGS-84, units=Degrees}' + '\n'
    sLine = 'wavelength units = Unknown '  + '\n'
    pFile.write(sLine)
    pFile.close()





