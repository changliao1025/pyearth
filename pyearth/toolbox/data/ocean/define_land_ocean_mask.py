#we will use gdal api for most operations
import os, sys
import numpy as np
import importlib.util
from pathlib import Path
from osgeo import ogr, osr, gdal
import cartopy.feature as cfeature

#https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.feature.NaturalEarthFeature.html

def create_land_ocean_vector_mask( sFilename_geojson_out ,
                           sResolution_coastal = '10m'):

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    # Create a land feature
    land_feature = cfeature.NaturalEarthFeature('physical', 'land', sResolution_coastal)
    land_geometries = list(land_feature.geometries())
    nPart = len(land_geometries)
    print(nPart)
    # Initialize pPolygon_land with the first geometry

    # Start the loop from the second geometry
    iFlag_first = 1
    iCount = 0
    sFolder = '/compyfs/liao313/04model/pyflowline/natural_earth'
    Path(sFolder).mkdir(parents=True, exist_ok=True)

    for land_geometry in land_geometries:
        pGeometry_mesh = ogr.CreateGeometryFromWkb(land_geometry.wkb)
        envelope = pGeometry_mesh.GetEnvelope()  # Returns a tuple (minX, maxX, minY, maxY)
        #check whether it is a multipolygon
        if pGeometry_mesh.GetGeometryName() == 'MULTIPOLYGON':
            nPart1 = pGeometry_mesh.GetGeometryCount()
            print(nPart1)
            for iPart in range(0, nPart1,1):
                pGeometry_part = pGeometry_mesh.GetGeometryRef(iPart)
                envelope = pGeometry_part.GetEnvelope()
                if envelope[2] < -60:
                    continue
                if iFlag_first == 1:
                    pPolygon_land = pGeometry_part
                    iFlag_first = 0
                else:
                    pPolygon_land = pPolygon_land.Union(pGeometry_part)

                sFilename_part = os.path.join(sFolder, f'land_geometry_{iCount}.geojson')
                if os.path.exists(sFilename_part):
                    os.remove(sFilename_part)

                pDataset_part = pDriver_geojson.CreateDataSource(sFilename_part)
                pLayer_part = pDataset_part.CreateLayer('part', geom_type=ogr.wkbPolygon)
                pFeature_part = ogr.Feature(pLayer_part.GetLayerDefn())
                pFeature_part.SetGeometry(pGeometry_part)
                pLayer_part.CreateFeature(pFeature_part)
                pFeature_part.Destroy()
                pDataset_part.Destroy()
                iCount = iCount + 1
        else:
            # Check if the geometry is in the region of Antarctica
            if envelope[2] < -60:
                # Skip this geometry
                #print(envelope)
                continue
            if iFlag_first == 1:
                pPolygon_land = ogr.CreateGeometryFromWkb(land_geometry.wkb)
                iFlag_first = 0
            else:
                pPolygon_land = pPolygon_land.Union(pGeometry_mesh)

            sFilename_part = os.path.join(sFolder, f'land_geometry_{iCount}.geojson')
            if os.path.exists(sFilename_part):
                os.remove(sFilename_part)

            pDataset_part = pDriver_geojson.CreateDataSource(sFilename_part)
            pLayer_part = pDataset_part.CreateLayer('part', geom_type=ogr.wkbPolygon)
            pFeature_part = ogr.Feature(pLayer_part.GetLayerDefn())
            pFeature_part.SetGeometry(pGeometry_mesh)
            pLayer_part.CreateFeature(pFeature_part)
            pFeature_part.Destroy()
            pDataset_part.Destroy()
            iCount = iCount + 1

    #save as output file
    if os.path.exists(sFilename_geojson_out):
        os.remove(sFilename_geojson_out)

    pDataset = pDriver_geojson.CreateDataSource(sFilename_geojson_out)
    pLayer = pDataset.CreateLayer('land_ocean_mask', geom_type=ogr.wkbPolygon)
    pFeature = ogr.Feature(pLayer.GetLayerDefn())
    pFeature.SetGeometry(pPolygon_land)
    pLayer.CreateFeature(pFeature)
    pFeature.Destroy()
    pDataset.Destroy()

    return

