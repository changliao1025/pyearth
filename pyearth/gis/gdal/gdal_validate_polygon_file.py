import os
from osgeo import ogr
def gdal_validate_polygon_file(sFilename_in):
    """
    Check if given polygon file is valid
    Args:
        sFilename_in (str): Filename of polygon file

    Returns:
        int: 1-valid, 0-invalid
    """

    iFlag_valid = 1

    # Check if file exists
    if not os.path.exists(sFilename_in):
        iFlag_valid = 0
    else:
        # Check if file is a polygon file
        try:
            pDataset = ogr.Open(sFilename_in)
            pLayer = pDataset.GetLayer()
            pFeature = pLayer.GetNextFeature()
            i = 0
            while pFeature is not None:
                #check if geometry is valid
                pGeometry = pFeature.GetGeometryRef()
                #check whether the geometry is a polygon
                if pGeometry.GetGeometryName() != 'POLYGON':
                    iFlag_valid = 0
                    break
                else:
                    #check whether the polygon is valid
                    if pGeometry.IsValid() == 0:
                        print("Polygon WKT:", pGeometry.ExportToWkt())
                        print(i)
                        iFlag_valid = 0
                        break
                    else:
                        pass
                i = i + 1
                pFeature = pLayer.GetNextFeature()
        except:
            iFlag_valid = 0

    return iFlag_valid
