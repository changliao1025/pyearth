import ogr, osr

def convert_xy_to_latlon(dX_in, dY_in, pSpatialRef_in):

    
    outputEPSG = 4326

    # create a geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(dX_in, dY_in)

    # create coordinate transformation
    #inSpatialRef = osr.SpatialReference()
    #inSpatialRef.ImportFromEPSG(inputEPSG)


    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outputEPSG)

    coordTransform = osr.CoordinateTransformation(pSpatialRef_in, outSpatialRef)

    # transform point
    point.Transform(coordTransform)

    

    dLongitude = point.GetX()
    dLatitude = point.GetY()

    return dLongitude, dLatitude