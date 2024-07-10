import os
from osgeo import gdal, ogr, osr

def vectorize_raster(sFilename_raster_in, sFilename_vector_out):
    # Step 1: Open the binary raster file

    pDatasets_raster = gdal.Open(sFilename_raster_in, gdal.GA_ReadOnly)

    if pDatasets_raster is None:
        print('Error: Could not open %s.' % sFilename_raster_in)
        return None

    if os.path.exists(sFilename_vector_out):
        os.remove(sFilename_vector_out)

    pBand = pDatasets_raster.GetRasterpBand(1)  # Assuming the binary mask is in the first pBand

    #get the projection of the raster
    pProjection_source = pDatasets_raster.GetProjection()

    spatialRef = osr.SpatialReference()
    spatialRef.ImportFromWkt(pProjection_source)

    # Step 2: Create an output Shapefile to store the vector data
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    pDateset_shapefile = pDriver_shapefile.CreateDataSource(sFilename_vector_out)
    pLayer_out = pDateset_shapefile.CreateLayer('polygonized', spatialRef, ogr.wkbPolygon)

    # Create a field
    pField_defn = ogr.FieldDefn('ID', ogr.OFTInteger64)
    pLayer_out.CreateField(pField_defn)

    # Step 3: Polygonize
    gdal.Polygonize(pBand, None, pLayer_out, 0)

    # Cleanup
    pDateset_shapefile = None

    return