import os
from osgeo import osr, ogr

def reproject_vector(sFilename_vector_in, sFilename_vector_out, pProjection_target):

    pDataset_vector = ogr.Open(sFilename_vector_in)    
    # Get the first layer in the shapefile
    pLayer_vector = pDataset_vector.GetLayer(0)   
    pSpatial_reference_vector = pLayer_vector.GetSpatialRef()
    pProjection_source = pSpatial_reference_vector.ExportToWkt()        

    pSpatial_reference_target = osr.SpatialReference()
    pSpatial_reference_target.ImportFromWkt(pProjection_target)

    pProjection_target = pSpatial_reference_target.ExportToWkt()

    #get the extension of polygon file  
    sExtension_vector = os.path.splitext(sFilename_vector_in)[1]
    #get the driver for the extension
    if sExtension_vector == '.geojson':
        pDriver_vector = ogr.GetDriverByName('GeoJSON')
    else:
        pDriver_vector = ogr.GetDriverByName(sExtension_vector)

    if (pProjection_source != pProjection_target):

        transform = osr.CoordinateTransformation(pSpatial_reference_vector, pSpatial_reference_target)

        pDataset_transform = pDriver_vector.CreateDataSource(sFilename_vector_out)
        #create a new shapefile layer
        pLayer_transform = pDataset_transform.CreateLayer('transform', pSpatial_reference_target, geom_type=ogr.wkbPolygon)
        #create a new feature
        pFeature_transform = ogr.Feature(pLayer_transform.GetLayerDefn())
        pFeature_clip = pLayer_vector.GetNextFeature()
        pPolygon = pFeature_clip.GetGeometryRef()
        pPolygon.Transform(transform)      
        #set the geometry
        pFeature_transform.SetGeometry(pPolygon)
        #add the feature to the layer
        pLayer_transform.CreateFeature(pFeature_transform)        
        pDataset_transform.FlushCache()
        #use the new shapefile to clip the raster
    else:
        print('The input and target vector files have the same spatial reference')    

    pSpatial_reference_vector = None
    pSpatial_reference_target = None
    return