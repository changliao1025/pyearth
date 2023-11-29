import os
from osgeo import ogr, gdal
def gdal_read_geojson_boundary(sFilename_boundary_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    iReturn_code = 1
    if os.path.isfile(sFilename_boundary_in):
        pass
    else:
        print('This mesh file does not exist: ', sFilename_boundary_in )
        iReturn_code = 0
        return iReturn_code

    
    pDriver_json = ogr.GetDriverByName('GeoJSON')    
    pDataset_mesh = pDriver_json.Open(sFilename_boundary_in, gdal.GA_ReadOnly)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_out = pLayer_mesh.GetSpatialRef()
    ldefn = pLayer_mesh.GetLayerDefn()   

    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()                     
        pGeometrytype_boundary = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_boundary == 'POLYGON'):       
            pBoundary_ogr = pGeometry_mesh  
        else:
            if(pGeometrytype_boundary == 'MULTIPOLYGON'):    
                nLine = pGeometry_mesh.GetGeometryCount()
                for i in range(nLine):
                    pBoundary_ogr = pGeometry_mesh.GetGeometryRef(i)
               
                pass
            else:
                pass
            pass   
            
            
    pBoundary_wkt = pBoundary_ogr.ExportToWkt()
    aExtent = pBoundary_ogr.GetEnvelope()
    min_x, max_x, min_y, max_y = aExtent
   
    return pBoundary_wkt, aExtent