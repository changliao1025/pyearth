from osgeo import ogr, osr, gdal

def clip_vector_by_shapefile(sFilename_vector_in, sFilename_polugon_in, sFilename_vector_out):
    # Open the input shapefile
    pDataset_source = ogr.Open(sFilename_vector_in)
    if pDataset_source is None:
        print("Error: Could not open the input shapefile.")
        return

    # Open the clip polygon
    pDataset_clip = ogr.Open(sFilename_polugon_in)
    if pDataset_clip is None:
        print("Error: Could not open the clip polygon.")
        return

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

    # Create a new layer in the output shapefile
    pLayer_clipped = pDataset_clipped.CreateLayer("clipped", geom_type=ogr.wkbPolygon)

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
            # Iterate over all features in the layer
        for feature in pLayer_clip:
            geometry = feature.GetGeometryRef()                
            if geometry is not None:
                geometry.Transform(transform)     
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
            new_feature = ogr.Feature(pLayer_clipped.GetLayerDefn())
            new_feature.SetGeometry(clipped_geometry)
            pLayer_clipped.CreateFeature(new_feature)
            new_feature = None

    # Close the shapefiles
    pDataset_source = None
    pDataset_clip = None
    pDataset_clipped = None

    print("Clipping completed.")

