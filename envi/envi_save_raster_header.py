def envi_save_raster_header(sFilename_in, aHeader_info_in):

    ofs = open(sFilename_in, 'w')
    
    for key, value in aHeader_info_in.items():
        sLine = key + ' = ' + value + '\n'
        ofs.write(sLine)

    ofs.close() 
  