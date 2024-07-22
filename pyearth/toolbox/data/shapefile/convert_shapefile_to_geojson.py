import os
from osgeo import ogr, gdal

def convert_shapefile_to_geojson(sFilename_shapefile_in, sFilename_geojson_out = None, sLayername_in = 'layer'):
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    #check if the pDriver_gpkg is available
    if pDriver_shapefile is None:
        print('shapefile driver not available.')
        print("GDAL version:", gdal.__version__)
        return


    pDataset_in = pDriver_shapefile.Open(sFilename_shapefile_in, 0)
    if pDataset_in is None:
        print(f"Failed to open file: {sFilename_shapefile_in}")
        return
    pLayer_in = pDataset_in.GetLayer()

    if sFilename_geojson_out is None:
        sFilename_geojson_out = sFilename_shapefile_in.replace('.shp', '.geojson')

    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_geojson_out)
    pLayer_out = pDataset_out.CopyLayer(pLayer_in, sLayername_in)
    pDataset_in = pDataset_out = None
    pLayer_in = pLayer_out = None

    return

