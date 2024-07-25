from osgeo import gdal
def gdal_check_raster_valid(sFilename_in, print_err=False):
    """Checks if `sFilename_in` is  valid file as tested from gdal's checksum function. Not an infalible check, but it is
    quick and should capture most cases.
    Adapted from https://lists.osgeo.org/pipermail/gdal-dev/2013-November/037520.html
    """
    sFilename_in = str(sFilename_in)
    try:
        ds = gdal.Open(sFilename_in)
        for i in range(ds.RasterCount):
            ds.GetRasterBand(i + 1).Checksum()
    except RuntimeError:
        valid = False
        if print_err:
            print('GDAL error: ', gdal.GetLastErrorMsg())
    else:
        valid = True
    return valid