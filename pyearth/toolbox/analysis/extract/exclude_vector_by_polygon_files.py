import os
from osgeo import ogr
from pyearth.toolbox.management.vector.reproject import reproject_vector
from pyearth.toolbox.management.vector.merge_features import merge_features

def exclude_vector_by_polygon_files(sFilename_vector_in, aFilename_polygon_in, sFilename_vector_out):
    # Open the input shapefile
    sExtension_source = os.path.splitext(sFilename_vector_in)[1]
    if sExtension_source == '.geojson':
        pDriver_source = ogr.GetDriverByName('GeoJSON')
    else:
        pDriver_source = ogr.GetDriverByName('ESRI Shapefile')
    pDataset_source = ogr.Open(sFilename_vector_in)
    if pDataset_source is None:
        print("Error: Could not open the input shapefile.")
        return

    sFolder_out = os.path.dirname(sFilename_vector_out)
    #delete the output file if it exists
    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)

    # Get the input layer and clip polygon layer
    pLayer_source = pDataset_source.GetLayer()
    #get the spatial reference
    pSpatial_reference_target = pLayer_source.GetSpatialRef()
    pProjection_target = pSpatial_reference_target.ExportToWkt()
    print(pProjection_target)

    # Create a new output
    pDataset_clipped = pDriver_source.CreateDataSource(sFilename_vector_out )
    if pDataset_clipped is None:
        print("Error: Could not create the output shapefile.")
        return

    #get the layer
    pLayer_in = pDataset_source.GetLayer()
    #obtain the geotype
    iGeomType = pLayer_in.GetGeomType()

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

    # Get the geometry of the clip polygon
    #get the extension of polygon file
    aFilename_polygon = list()
    for sFilename_polygon_in in aFilename_polygon_in:
        # Open the clip polygon
        pDataset_clip = ogr.Open(sFilename_polygon_in)
        if pDataset_clip is None:
            print("Error: Could not open the clip polygon.")
            return
        sExtension_clip = os.path.splitext(sFilename_vector_in)[1]
        #get the driver for the extension
        if sExtension_clip == '.geojson':
            pDriver_clip = ogr.GetDriverByName('GeoJSON')
        else:
            pDriver_clip = ogr.GetDriverByName('ESRI Shapefile')

        if pDriver_clip is None:
            print('The driver is not available')
            return

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
            sFilename_clip_out = sFolder + '/' + sName_no_extension + '_transformed' + sExtension_clip
            reproject_vector(sFilename_polygon_in, sFilename_clip_out, pProjection_target)
            pDataset_clip = ogr.Open(sFilename_clip_out)
            pLayer_clip = pDataset_clip.GetLayer(0)
            aFilename_polygon.append(sFilename_clip_out)
        else:
            aFilename_polygon.append(sFilename_polygon_in)

    #merge the polygon files into one file
    if len(aFilename_polygon) > 1:
        pGeometry_merge = ogr.Geometry(ogr.wkbPolygon)
        sFilenane_merge = os.path.join(sFolder_out, 'polygon_merged' + sExtension_source)
        if os.path.exists(sFilenane_merge):
            os.remove(sFilenane_merge)
        pDataset_merge = pDriver_clip.CreateDataSource(sFilenane_merge )
        if pDataset_merge is None:
            print("Error: Could not create the output shapefile.")
            return
        pLayer_merge = pDataset_merge.CreateLayer('layer', pSpatial_reference_target, geom_type=ogr.wkbMultiPolygon)
        iFlag_first = 1

        for sFilename_polygon_in in aFilename_polygon:
            pDataset_clip = ogr.Open(sFilename_polygon_in)
            pLayer_clip = pDataset_clip.GetLayer()
            for pFeature_clip in pLayer_clip:
                pLayer_merge.CreateFeature(pFeature_clip)

            pGeometry = pFeature_clip.GetGeometryRef()
            if iFlag_first == 1:
                if pGeometry is not None:
                    # Union the geometry of each feature with the merged polygon
                    pGeometry_merge = pGeometry_merge.Union(pGeometry)

                iFlag_first = 0
            else:
                if pGeometry is not None:
                    pGeometry_merge = pGeometry_merge.Union(pGeometry)
    else:
        # Open the clip polygon
        pDataset_clip = ogr.Open(aFilename_polygon[0])
        if pDataset_clip is None:
            print("Error: Could not open the clip polygon.")
            return
        pLayer_clip = pDataset_clip.GetLayer()
        # Get the geometry of the clip polygon using next feature
        pFeature_clip = pLayer_clip.GetNextFeature()
        pGeometry_merge = pFeature_clip.GetGeometryRef()

    #use the merged polygon as the clip polygon
    pPolygon_clip = pGeometry_merge
    #print its envelope
    aEnvelope = pPolygon_clip.GetEnvelope()
    print('Envelope:', aEnvelope)

    # Apply the clipping operation to each pFeature in the input shapefile
    for pFeature in pLayer_source:
        # Get the geometry of the input pFeature
        pGeometry_source = pFeature.GetGeometryRef()
        #check pFeature that is entirely within the clip polygon
        if pPolygon_clip.Contains(pGeometry_source):
            #if it is within, we definitely do not want to keep it
            pass
        else:
            # Perform the clipping operation
            iFlag_ok = 0
            for i in range(pGeometry_source.GetPointCount()):
                aPoint = pGeometry_source.GetPoint(i)
                pPoint = ogr.Geometry(ogr.wkbPoint)
                pPoint.AddPoint(*aPoint)
                iFlag_within = pPoint.Within(pPolygon_clip)
                if iFlag_within:
                    iFlag_ok = 1
                    break
            if iFlag_ok == 1:
                if pGeometry_source.Intersects(pPolygon_clip):
                    pGeometry_clip = pGeometry_source.Difference(pPolygon_clip)
                    #Create a new pFeature in the output layer with the clipped geometry
                    if pGeometry_clip is not None and not pGeometry_clip.IsEmpty():
                        iGeomType_intersect = pGeometry_clip.GetGeometryType()
                        if iGeomType_intersect == iGeomType:
                            #we only want to keep the clipped geometry that is outside the clip polygon
                            nPart = pGeometry_clip.GetGeometryCount()
                            #print('Number of parts:', nPart)
                            pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                            pFeature_clipped.SetGeometry(pGeometry_clip)
                            #print pGeometry_clip to wkt
                            print(pGeometry_clip.ExportToWkt())
                            for i in range(0, pFeature.GetFieldCount()):
                                pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))
                            pLayer_clipped.CreateFeature(pFeature_clipped)
                        else:
                            sGeomType = ogr.GeometryTypeToName(iGeomType_intersect)
                            print('Intersect type:', sGeomType)
                            #it is possible that the clipped geometry is a multi linesting
                            if iGeomType_intersect == ogr.wkbMultiLineString:
                                nPart = pGeometry_clip.GetGeometryCount()
                                print('Number of parts here:', nPart)
                                for iPart in range(nPart):
                                    pGeometry_part = pGeometry_clip.GetGeometryRef(iPart)
                                    pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                                    pFeature_clipped.SetGeometry(pGeometry_part)
                                    for i in range(0, pFeature.GetFieldCount()):
                                        pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))
                                    pLayer_clipped.CreateFeature(pFeature_clipped)
                            continue
            else:
                #we actually want to keep the pFeature because it is entirely outside the clip polygon
                pFeature_clipped = ogr.Feature(pLayer_clipped.GetLayerDefn())
                pFeature_clipped.SetGeometry(pGeometry_source)
                # Copy attributes
                for i in range(0, pFeature.GetFieldCount()):
                    pFeature_clipped.SetField(pFeature.GetFieldDefnRef(i).GetNameRef(), pFeature.GetField(i))
                pLayer_clipped.CreateFeature(pFeature_clipped)
                continue

    # Close the shapefiles
    pDataset_source = None
    pDataset_clip = None
    pDataset_clipped = None
    pSpatial_reference_clip = None
    pSpatial_reference_target = None

    print("Clipping completed.")