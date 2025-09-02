import os
import numpy as np
from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter, meter_to_degree



def find_nearest_resolution(dResolution_meter_in, dLatitude_in):

    #convert the resolution to degree
    dResolution_degree = meter_to_degree(dResolution_meter_in, dLatitude_in)
    #derive the resolution in between 1 , 1/2, 1/4, 1/8, 1/16, 1/32, 1/64
    nResolution = 10
    aResolution_degree = np.zeros((nResolution), dtype=float)
    for i in range(0, nResolution):
        aResolution_degree[i] = 1.0 / (2**i)

    #find the nearest resolution that is larger than the input resolution
    index = 0
    for i in range(0, nResolution-1):
        if aResolution_degree[i] >= dResolution_degree and aResolution_degree[i+1] < dResolution_degree:
            index = i
            break

    dResolution_degree_out = aResolution_degree[index]

    #double check
    dRes = degree_to_meter(dResolution_degree_out, dLatitude_in)

    return dResolution_degree_out

def create_box_from_longitude_latitude(dLongitude_in, dLatitude_in,
                                       dResolution_x_in, dResolution_y_in, sFilename_output_in=None):


    if sFilename_output_in is not None and os.path.isfile(sFilename_output_in):
        print(sFilename_output_in + ' already exists')
        os.remove(sFilename_output_in)



    #dLon_min = dLongitude_in - dResolution_x_in/2
    #dLon_max = dLongitude_in + dResolution_x_in/2
    #dLat_min = dLatitude_in - dResolution_y_in/2
    #dLat_max = dLatitude_in + dResolution_y_in/2

    dLon_min = dLongitude_in
    dLat_max = dLatitude_in


    #determine structure boundary
    #use the global as the
    #the resolution must be 0.5, 0.25 etc
    #maybe adding a checkpoint here to confirm it

    nleft  = np.floor(  (dLon_min - (-180)) /(dResolution_x_in)  )
    ntop  = np.floor(  (90 - dLat_max) /(dResolution_y_in)  )

    #should be 1
    aBox_out = np.full( (4), -9999, dtype=float )
    aCoordinates_out = np.full( (4,2), -9999, dtype=float )


    #left
    aBox_out[0] = -180 + (nleft ) * dResolution_x_in #+ 0.5 * dResolution_x_in
    #top
    aBox_out[3] = 90 - ( ntop )  * dResolution_y_in  #- 0.5 * dResolution_y_in
    #right
    aBox_out[1] = aBox_out[0] + dResolution_x_in
    #bot
    aBox_out[2] = aBox_out[3] - dResolution_y_in

    #contercolckwise
    #upper left
    aCoordinates_out[0,0] = aBox_out[0]
    aCoordinates_out[0,1] = aBox_out[3]

    #lower left
    aCoordinates_out[1,0] = aBox_out[0]
    aCoordinates_out[1,1] = aBox_out[2]

    #lower right
    aCoordinates_out[2,0] = aBox_out[1]
    aCoordinates_out[2,1] = aBox_out[2]

    #upper right
    aCoordinates_out[3,0] = aBox_out[1]
    aCoordinates_out[3,1] = aBox_out[3]

    dArea = calculate_polygon_area(aCoordinates_out[:,0], aCoordinates_out[:,1])

    #save as a geojson file
    if sFilename_output_in is not None:

        pDriver_geojson = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
        pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
        # Add one attribute
        pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
        pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
        pArea_field.SetWidth(20)
        pArea_field.SetPrecision(2)
        pLayer.CreateField(pArea_field)
        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in aCoordinates_out:
            ring.AddPoint(x, y)

        ring.AddPoint(aCoordinates_out[0,0], aCoordinates_out[0,1])

        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        pFeature.SetGeometry(pPolygon)
        pFeature.SetField("cellid", 1)
        pFeature.SetField("longitude", dLongitude_in )
        pFeature.SetField("latitude", dLatitude_in )
        pFeature.SetField("area", dArea )
        pLayer.CreateFeature(pFeature)

        pDataset = pLayer = pFeature  = None

    return aBox_out, aCoordinates_out