from osgeo import gdal, osr, ogr
def gdal_get_vector_extent(sFilename_in):
    #get geojson driver

    pDataset_in = ogr.Open(sFilename_in)
    # Get the first layer in the file
    pLayer_in = pDataset_in.GetLayer(0)
    extent = pLayer_in.GetExtent()
    #get min_x, max_x, min_y, max_y
    min_x = extent[0]
    max_x = extent[1]
    min_y = extent[2]
    max_y = extent[3]
    return extent