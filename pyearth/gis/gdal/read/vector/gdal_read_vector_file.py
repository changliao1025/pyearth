import os
from osgeo import gdal, ogr
import numpy as np


def gdal_read_shapefile(sFilename_in, iFlag_metadata_only=0):
    """
    Read a ESRI Shapefile.
    Be careful that GDAL Python blinding has issue with some objects

    Args:
        sFilename_in (string): The filename

    Returns:
        tuple: aFeature, pSpatial_reference, lFeatureCount
    """

    if os.path.exists(sFilename_in):
        pass
    else:
        print("The file does not exist!")
        return

    sDriverName = "ESRI Shapefile"
    pDriver = ogr.GetDriverByName(sDriverName)

    if pDriver is None:
        print("%s pDriver not available.\n" % sDriverName)
    else:
        print("%s pDriver IS available.\n" % sDriverName)

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)

    # Check to see if shapefile is found.
    if pDataset is None:
        print("Could not open %s" % (sFilename_in))
    else:
        pLayer = pDataset.GetLayer()
        lFeatureCount = pLayer.GetFeatureCount()
        print(
            "Number of features in %s: %d"
            % (os.path.basename(sFilename_in), lFeatureCount)
        )

        pSpatial_reference = pLayer.GetSpatialRef()
        iFlag_first = 1

        for feature in pLayer:
            pGeometry = feature.GetGeometryRef()
            pEnvelope = pGeometry.GetEnvelope()

            if iFlag_first == 1:
                left_min = pEnvelope[0]
                right_max = pEnvelope[1]
                bot_min = pEnvelope[2]
                top_max = pEnvelope[3]
                iFlag_first = 0
            else:
                left_min = np.min([left_min, pEnvelope[0]])
                right_max = np.max([right_max, pEnvelope[1]])
                bot_min = np.min([bot_min, pEnvelope[2]])
                top_max = np.max([top_max, pEnvelope[3]])

        if iFlag_metadata_only == 1:
            return (
                pSpatial_reference,
                lFeatureCount,
                left_min,
                right_max,
                bot_min,
                top_max,
            )
        else:
            aFeature = []
            for pFeature in pLayer:
                pGeometry = pFeature.GetGeometryRef()
                aFeature.append(pGeometry)
                pass

            return (
                aFeature,
                pSpatial_reference,
                lFeatureCount,
                left_min,
                right_max,
                bot_min,
                top_max,
            )
