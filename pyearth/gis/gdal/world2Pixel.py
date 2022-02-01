
def world2Pixel(pGeoMatrix_in, dx_in, dy_in):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
   
    Args:
        geoMatrix (gdal): The geotransform matrix
        dx_in ([type]): The X coodinate of a point
        dy_in ([type]): The Y coodinate of a point

    Returns:
        Tuple: pixel, line
    """
    
    ulX = pGeoMatrix_in[0]
    ulY = pGeoMatrix_in[3]
    xDist = pGeoMatrix_in[1]
    yDist = pGeoMatrix_in[5]
    rtnX = pGeoMatrix_in[2]
    rtnY = pGeoMatrix_in[4]
    pixel = int((dx_in - ulX) / xDist)
    line = int((ulY - dy_in) / xDist)
    return (pixel, line)