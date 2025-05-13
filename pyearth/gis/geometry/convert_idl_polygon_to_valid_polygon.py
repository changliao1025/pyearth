import numpy as np
from osgeo import ogr
from pyearth.system.define_global_variables import *
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
def convert_idl_polygon_to_valid_polygon(pGeometry_in):

    #get the coordinates first
    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
    dLongitude_min= np.min(aCoords_gcs[:,0])
    dLongitude_max = np.max(aCoords_gcs[:,0])
    if dLongitude_max - dLongitude_min < 180:
        if pGeometry_in.IsValid():
            return pGeometry_in
        else:
            return None
    else:
        nPoints = len(aCoords_gcs)
        for i in range(nPoints):
            dLongitude = aCoords_gcs[i][0]
            if dLongitude < -150 :
                dLongitude = dLongitude + 360.0
                pass
            aCoords_gcs[i][0] = dLongitude
            pass
        pGeometry_out = ogr.Geometry(ogr.wkbPolygon)
        pRing = ogr.Geometry(ogr.wkbLinearRing)
        for aCoord in aCoords_gcs:
            pRing.AddPoint(aCoord[0], aCoord[1])
            pass
        pRing.CloseRings()
        pGeometry_out.AddGeometry(pRing)
        pGeometry_out.AssignSpatialReference(pGeometry_in.GetSpatialReference())
        #check valid
        if pGeometry_out.IsValid() == False:
            pGeometry_out = None
            pass
        return pGeometry_out

if __name__ == '__main__':
    #check the validity of a polygon if it cross the international date line
    wkt = 'POLYGON ((-179.997342114428 66.2769360301573,179.98848859037 66.2814388694679,179.984323361886 66.281448567623,179.976318838375 66.27728074511330,179.976930616892 66.275689036998,179.98551001155 66.2729876032834,179.990204009696 66.2729766264868,-179.997342114428 66.2769360301573))'
    pGeometry_in = ogr.CreateGeometryFromWkt(wkt)
    if pGeometry_in is None:
        print('failed to create geometry')
        pass
    else:
        if pGeometry_in.IsValid():
            print('The input geometry is valid')
        else:
            print('The input geometry is invalid')

            pGeometry_out = convert_idl_polygon_to_valid_polygon(pGeometry_in)
            if pGeometry_out is None:
                print('The output geometry is invalid')
            else:
                if pGeometry_out.IsValid():
                    print('The output geometry is valid')
                else:
                    print('The output geometry is invalid')
                    pass

