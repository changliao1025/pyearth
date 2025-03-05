import os
from osgeo import ogr
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.toolbox.management.vector.merge_features import merge_features

def clip_vector_by_polygon_file(sFilename_vector_in, sFilename_polygon_in, sFilename_vector_out ):
    # Open the input shapefile
    #check file existence
    if os.path.exists(sFilename_vector_in) == False:
        print(f"Error: file {sFilename_vector_in} does not exist!")
        return

    pDataset_source = ogr.Open(sFilename_vector_in)
    if pDataset_source is None:
        print("Error: Could not open the input shapefile.")
        return

    # Open the clip polygon
    pDataset_clip = ogr.Open(sFilename_polygon_in)
    if pDataset_clip is None:
        print("Error: Could not open the clip polygon.")
        return

    #get the extension of polygon file
    sExtension_vector = os.path.splitext(sFilename_polygon_in)[1]
    #get the driver for the extension
    if sExtension_vector == '.geojson':
        pDriver_vector = ogr.GetDriverByName('GeoJSON')
    else:
        pDriver_vector = ogr.GetDriverByName('ESRI Shapefile')

    if pDriver_vector is None:
        print('The driver is not available')
        return

    #delete the output file if it exists
    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)

    # Get the input layer and clip polygon layer
    pLayer_source = pDataset_source.GetLayer()
    #get the spatial reference
    pSpatial_reference_target = pLayer_source.GetSpatialRef()
    pProjection_target = pSpatial_reference_target.ExportToWkt()
    print(pProjection_target)

    # Get the geometry of the clip polygon
    sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]
    pLayer_clip = pDataset_clip.GetLayer()
    pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
    pProjection_clip = pSpatial_reference_clip.ExportToWkt()
    print(pProjection_clip)
    #check the number the clip polygon
    nPolygon = pLayer_clip.GetFeatureCount()
    if nPolygon == 0:
        print('The polygon file does not contain any polygon!')
        return
    else:
        if nPolygon > 1:
            pDataset_clip = None
            pLayer_clip = None
            print('The polygon contains more than one polygon, the program will attempt to merge them as one!')
            #obtain the file extension
            sFilename_clip_new = sFilename_polygon_in.replace(sExtension_clip, '_merged' + sExtension_clip)
            merge_features(sFilename_polygon_in, sFilename_clip_new)
            sFilename_polygon_in = sFilename_clip_new
            #open the new file
            pDataset_clip = ogr.Open(sFilename_polygon_in)
            pLayer_clip = pDataset_clip.GetLayer(0)
            # Get the spatial reference of the layer
            pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
            pProjection_clip = pSpatial_reference_clip.ExportToWkt()
        else:
            pass

    if( pProjection_target != pProjection_clip):
        pDataset_clip = None
        pLayer_clip = None
        sFolder = os.path.dirname(sFilename_polygon_in)
        #get the name of the shapefile
        sName = os.path.basename(sFilename_polygon_in)
        #get the name of the shapefile without extension
        sName_no_extension = os.path.splitext(sName)[0]
        sFilename_clip_out = sFolder + '/' + sName_no_extension + '_transformed' + sExtension_vector
        reproject_vector(sFilename_polygon_in, sFilename_clip_out, pProjection_target)
        pDataset_clip = ogr.Open(sFilename_clip_out)
        pLayer_clip = pDataset_clip.GetLayer(0)

    #the the extent rectangle of the clip polygon
    pEnvelope = pLayer_clip.GetExtent()
    minx, maxx, miny,  maxy = pEnvelope
    pLayer_source.SetSpatialFilterRect(minx, miny, maxx, maxy)

    pLayer_clip.ResetReading()
    pFeature_clip = pLayer_clip.GetNextFeature()
    pPolygon_clip = pFeature_clip.GetGeometryRef()

    # Create a new output
    pDataset_clipped = pDriver_vector.CreateDataSource(sFilename_vector_out )
    if pDataset_clipped is None:
        print("Error: Could not create the output shapefile.")
        return

    #get the layer
    #obtain the geotype
    #get the first layer
    pLayer_in = pDataset_source.GetLayer()
    #get the first feature
    pFeature_in = pLayer_in.GetNextFeature()
    #get the first geometry
    pGeometry_in = pFeature_in.GetGeometryRef()
    #get the geometry type
    iGeomType = pGeometry_in.GetGeometryType()
    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString :
            #create the layer
            pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPolygon)
            else:
                if iGeomType == ogr.wkbLineString25D:
                    pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                    iGeomType = ogr.wkbLineString
                else:
                    if iGeomType == ogr.wkbMultiPolygon:
                        pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPolygon)
                        iGeomType = ogr.wkbPolygon
                    else:
                        if iGeomType == ogr.wkbMultiLineString:
                            pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                            iGeomType = ogr.wkbLineString
                        else:
                            if iGeomType == ogr.wkbMultiLineString25D:
                                pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                                iGeomType = ogr.wkbLineString
                            else:
                                sGeomType = ogr.GeometryTypeToName(iGeomType)
                                print('Geometry type not supported:', sGeomType)
                                return
            pass
        pass


    # Get the feature definition of the source layer
    pFeatureDefn_source = pLayer_source.GetLayerDefn()

    # Create the fields in the clipped layer
    for i in range(pFeatureDefn_source.GetFieldCount()):
        pFieldDefn_source = pFeatureDefn_source.GetFieldDefn(i)
        pLayer_clipped.CreateField(pFieldDefn_source)

    # Apply the clipping operation to each pFeature in the input shapefile
    lID =1
    for pFeature in pLayer_source:
        # Get the geometry of the input pFeature
        pGeometry_source = pFeature.GetGeometryRef()
        #check pFeature that is entirely within the clip polygon
        if pPolygon_clip.Contains(pGeometry_source):
            pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
            pFeature_clipped.SetGeometry(pGeometry_source)
            # Copy attributes
            for i in range(0, pFeature.GetFieldCount()):
                pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))

            pLayer_clipped.CreateFeature(pFeature_clipped)
        else:
            # Perform the clipping operation
            if pGeometry_source.Intersects(pPolygon_clip) == True:
                pGeometry_intersect = pGeometry_source.Intersection(pPolygon_clip)
                if pGeometry_intersect is None or pGeometry_intersect.IsEmpty():
                    print("Intersection resulted in an empty or invalid geometry")
                    continue
                else:
                    iGeomType_intersect = pGeometry_intersect.GetGeometryType()
                    print(f"Intersection geometry type: {iGeomType_intersect}")

                    # Create a new pFeature in the output layer with the clipped geometry
                    #we only want to keep the clipped geometry that is within the clip polygon
                    pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                    pFeature_clipped.SetGeometry(pGeometry_intersect)
                    for i in range(0, pFeature.GetFieldCount()):
                        pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))
                    pLayer_clipped.CreateFeature(pFeature_clipped)



    # Close the shapefiles
    pDataset_source = None
    pDataset_clip = None
    pDataset_clipped = None
    pSpatial_reference_clip = None
    pSpatial_reference_target = None

    print("Clipping completed.")

#multipe file to clip
def clip_vector_by_polygon_files(sFilename_vector_in, aFilename_polygon_in, sFilename_vector_out):
    # Open the input shapefile
    pDataset_source = ogr.Open(sFilename_vector_in)
    if pDataset_source is None:
        print("Error: Could not open the input shapefile.")
        return

    #delete the output file if it exists
    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)

    # Open the clip polygon
    pDataset_clip = ogr.Open(sFilename_polygon_in)
    if pDataset_clip is None:
        print("Error: Could not open the clip polygon.")
        return

    # Get the input layer and clip polygon layer
    pLayer_source = pDataset_source.GetLayer()
    #get the spatial reference
    pSpatial_reference_target = pLayer_source.GetSpatialRef()
    pProjection_target = pSpatial_reference_target.ExportToWkt()
    print(pProjection_target)

    # Create a new output
    pDataset_clipped = pDriver_vector.CreateDataSource(sFilename_vector_out )
    if pDataset_clipped is None:
        print("Error: Could not create the output shapefile.")
        return

    #get the layer
    #obtain the geotype
    #get the first layer
    pLayer_in = pDataset_source.GetLayer()
    #get the first feature
    pFeature_in = pLayer_in.GetNextFeature()
    #get the first geometry
    pGeometry_in = pFeature_in.GetGeometryRef()
    #get the geometry type
    iGeomType = pGeometry_in.GetGeometryType()
    if iGeomType == ogr.wkbPoint:
        #create the layer
        pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPoint)
    else:
        if iGeomType == ogr.wkbLineString :
            #create the layer
            pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
        else:
            if iGeomType == ogr.wkbPolygon:
                #create the layer
                pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPolygon)
            else:
                if iGeomType == ogr.wkbLineString25D:
                    pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                    iGeomType = ogr.wkbLineString
                else:
                    if iGeomType == ogr.wkbMultiPolygon:
                        pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbPolygon)
                        iGeomType = ogr.wkbPolygon
                    else:
                        if iGeomType == ogr.wkbMultiLineString:
                            pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                            iGeomType = ogr.wkbLineString
                        else:
                            if iGeomType == ogr.wkbMultiLineString25D:
                                pLayer_clipped = pDataset_clipped.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbLineString)
                                iGeomType = ogr.wkbLineString
                            else:
                                sGeomType = ogr.GeometryTypeToName(iGeomType)
                                print('Geometry type not supported:', sGeomType)
                                return
            pass
        pass

    # Get the feature definition of the source layer
    pFeatureDefn_source = pLayer_source.GetLayerDefn()

    # Create the fields in the clipped layer
    for i in range(pFeatureDefn_source.GetFieldCount()):
        pFieldDefn_source = pFeatureDefn_source.GetFieldDefn(i)
        pLayer_clipped.CreateField(pFieldDefn_source)

    #get the extension of polygon file
    for sFilename_polygon_in in aFilename_polygon_in:
        sExtension_vector = os.path.splitext(sFilename_vector_in)[1]
        #get the driver for the extension
        if sExtension_vector == '.geojson':
            pDriver_vector = ogr.GetDriverByName('GeoJSON')
        else:
            pDriver_vector = ogr.GetDriverByName('ESRI Shapefile')

        if pDriver_vector is None:
            print('The driver is not available')
            return


        # Get the geometry of the clip polygon
        sExtension_clip = os.path.splitext(sFilename_polygon_in)[1]
        pLayer_clip = pDataset_clip.GetLayer()
        pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
        pProjection_clip = pSpatial_reference_clip.ExportToWkt()
        print(pProjection_clip)
        #check the number the clip polygon
        nPolygon = pLayer_clip.GetFeatureCount()
        if nPolygon == 0:
            print('The polygon file does not contain any polygon!')
            return
        else:
            if nPolygon > 1:
                pDataset_clip = None
                pLayer_clip = None
                print('The polygon contains more than one polygon, the program will attempt to merge them as one!')
                #obtain the file extension
                sFilename_clip_new = sFilename_polygon_in.replace(sExtension_clip, '_merged' + sExtension_clip)
                merge_features(sFilename_polygon_in, sFilename_clip_new)
                sFilename_polygon_in = sFilename_clip_new
                #open the new file
                pDataset_clip = ogr.Open(sFilename_polygon_in)
                pLayer_clip = pDataset_clip.GetLayer(0)
                # Get the spatial reference of the layer
                pSpatial_reference_clip = pLayer_clip.GetSpatialRef()
                pProjection_clip = pSpatial_reference_clip.ExportToWkt()
            else:
                pass

        if( pProjection_target != pProjection_clip):
            pDataset_clip = None
            pLayer_clip = None
            sFolder = os.path.dirname(sFilename_polygon_in)
            #get the name of the shapefile
            sName = os.path.basename(sFilename_polygon_in)
            #get the name of the shapefile without extension
            sName_no_extension = os.path.splitext(sName)[0]
            sFilename_clip_out = sFolder + '/' + sName_no_extension + '_transformed' + sExtension_vector
            reproject_vector(sFilename_polygon_in, sFilename_clip_out, pProjection_target)
            pDataset_clip = ogr.Open(sFilename_clip_out)
            pLayer_clip = pDataset_clip.GetLayer(0)

        #the the extent rectangle of the clip polygon
        pEnvelope = pLayer_clip.GetExtent()
        minx, maxx, miny,  maxy = pEnvelope
        pLayer_source.SetSpatialFilterRect(minx, miny, maxx, maxy)

        pLayer_clip.ResetReading()
        pFeature_clip = pLayer_clip.GetNextFeature()
        pPolygon_clip = pFeature_clip.GetGeometryRef()

        # Apply the clipping operation to each pFeature in the input shapefile
        for pFeature in pLayer_source:
            # Get the geometry of the input pFeature
            pGeometry_source = pFeature.GetGeometryRef()
            #check pFeature that is entirely within the clip polygon
            if pPolygon_clip.Contains(pGeometry_source):
                pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                pFeature_clipped.SetGeometry(pGeometry_source)
                # Copy attributes
                for i in range(0, pFeature.GetFieldCount()):
                    pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))

                pLayer_clipped.CreateFeature(pFeature_clipped)
            else:
                # Perform the clipping operation
                # Perform the clipping operation
                pGeometry_intersect = pGeometry_source.Intersection(pPolygon_clip)
                if pGeometry_intersect is None or pGeometry_intersect.IsEmpty():
                    print("Intersection resulted in an empty or invalid geometry")
                    continue
                else:
                    iGeomType_intersect = pGeometry_intersect.GetGeometryType()
                    print(f"Intersection geometry type: {iGeomType_intersect}")

                    # Create a new pFeature in the output layer with the clipped geometry
                    #we only want to keep the clipped geometry that is within the clip polygon
                    pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                    pFeature_clipped.SetGeometry(pGeometry_intersect)
                    for i in range(0, pFeature.GetFieldCount()):
                        pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))
                    pLayer_clipped.CreateFeature(pFeature_clipped)

    # Close the shapefiles
    pDataset_source = None
    pDataset_clip = None
    pDataset_clipped = None
    pSpatial_reference_clip = None
    pSpatial_reference_target = None

    print("Clipping completed.")