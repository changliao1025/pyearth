import os
from osgeo import ogr, osr, gdal

def clip_vector_by_shapefile(sFilename_vector_in, sFilename_polygon_in, sFilename_vector_out):
    # Open the input shapefile
    pDataset_source = ogr.Open(sFilename_vector_in)
    if pDataset_source is None:
        print("Error: Could not open the input shapefile.")
        return

    # Open the clip polygon
    pDataset_clip = ogr.Open(sFilename_polygon_in)
    if pDataset_clip is None:
        print("Error: Could not open the clip polygon.")
        return
    
    #delete the output shapefile file if it exists
    driver = ogr.GetDriverByName("ESRI Shapefile")

    if os.path.exists(sFilename_vector_out):
        driver.DeleteDataSource(sFilename_vector_out)
   

    # Get the input layer and clip polygon layer
    pLayer_source = pDataset_source.GetLayer()
    
    #get the spatial reference
    pSpatial_reference_source = pLayer_source.GetSpatialRef()


    # Get the geometry of the clip polygon
    pLayer_clip = pDataset_clip.GetLayer()
    
    pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
    comparison = pSpatial_reference_source.IsSame(pSpatial_reference_clip)
    
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_source, pSpatial_reference_clip)
        pass
    else:
        iFlag_transform = 0
        pass

    # Create a new output shapefile
    driver = ogr.GetDriverByName("ESRI Shapefile")
    pDataset_clipped = driver.CreateDataSource(sFilename_vector_out)
    if pDataset_clipped is None:
        print("Error: Could not create the output shapefile.")
        return

    #get the layer
    pLayer_in = pDataset_source.GetLayer()
    #obtain the geotype 
    iGeomType = pLayer_in.GetGeomType()
    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_source, geom_type=ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString:
            #create the layer
            pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_source, geom_type=ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_source, geom_type=ogr.wkbPolygon)
            else:
                print('Geometry type not supported')
                return
            pass
        pass
    # Create a new layer in the output shapefile
    

    # Define the clipping operation
  
    #check the number the clip polygon
    iCount = pLayer_clip.GetFeatureCount()
    if iCount ==1:        
        pFeature = pLayer_clip.GetNextFeature()
        pPolygon_clip = pFeature.GetGeometryRef()
        pass
    else:
        #merge the clip polygon
        pPolygon_clip = ogr.Geometry(ogr.wkbPolygon)

        if iFlag_transform == 1:
            # Iterate over all features in the layer
            for feature in pLayer_clip:
                geometry = feature.GetGeometryRef()                
                if geometry is not None:
                    geometry.Transform(transform)     
                    # Union the geometry of each feature with the merged polygon
                    pPolygon_clip = pPolygon_clip.Union(geometry)     
        else:
            # Iterate over all features in the layer
            for feature in pLayer_clip:
                geometry = feature.GetGeometryRef()
                if geometry is not None:
                    # Union the geometry of each feature with the merged polygon
                    pPolygon_clip = pPolygon_clip.Union(geometry)        
  

    # Apply the clipping operation to each feature in the input shapefile
    for feature in pLayer_source:
        # Get the geometry of the input feature
        input_geometry = feature.GetGeometryRef()

        # Perform the clipping operation
        clipped_geometry = input_geometry.Intersection(pPolygon_clip)

        # Create a new feature in the output layer with the clipped geometry
        if clipped_geometry is not None and not clipped_geometry.IsEmpty():

            #we only want to keep the clipped geometry that is within the clip polygon
            
            new_feature = ogr.Feature(pLayer_clipped.GetLayerDefn())
            new_feature.SetGeometry(clipped_geometry)
            pLayer_clipped.CreateFeature(new_feature)
            new_feature = None
            print('intersected')

        else:

            #check feature that is entirely within the clip polygon
            if pPolygon_clip.Contains(input_geometry):
                new_feature = ogr.Feature(pLayer_clipped.GetLayerDefn())
                new_feature.SetGeometry(input_geometry)
                pLayer_clipped.CreateFeature(new_feature)
                new_feature = None
                print('within')
        


    # Close the shapefiles
    pDataset_source = None
    pDataset_clip = None
    pDataset_clipped = None

    print("Clipping completed.")
