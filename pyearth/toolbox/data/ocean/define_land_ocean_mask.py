#we will use gdal api for most operations
import os, sys
import numpy as np
import importlib.util
from pathlib import Path
from osgeo import ogr, osr, gdal
import cartopy.feature as cfeature
from pyearth.gis.gdal.read.vector.get_supported_formats import get_format_from_extension

#https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.feature.NaturalEarthFeature.html



def create_land_ocean_vector_mask(sFilename_out ,
                           sResolution_coastal = '10m',
                           sWorkspace_out = None,
                           sFormat = None):
    """
    Create a land-ocean vector mask from Natural Earth data.

    Args:
        sFilename_out: Output filename (format determined by extension if sFormat not specified)
        sResolution_coastal: Resolution for coastal data ('10m', '50m', '110m')
        sWorkspace_out: Optional workspace folder for intermediate files
        sFormat: Output format ('GeoJSON', 'ESRI Shapefile', 'GPKG', etc.)
    """

    # Auto-detect format from output filename if not specified
    if sFormat is None:
        sFormat = get_format_from_extension(sFilename_out)
        print(f'Auto-detected format: {sFormat}')

    # Get the appropriate driver
    pDriver_out = ogr.GetDriverByName(sFormat)
    if pDriver_out is None:
        print(f'{sFormat} driver not available.')
        return

    # Keep GeoJSON driver for intermediate files
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    # Create a land feature
    land_feature = cfeature.NaturalEarthFeature('physical', 'land', sResolution_coastal)
    land_geometries = list(land_feature.geometries())
    nPart = len(land_geometries)
    print(nPart)

    # Create output dataset first
    if os.path.exists(sFilename_out):
        # For shapefiles, we need to delete all associated files
        if sFormat == 'ESRI Shapefile':
            # Get the base name without extension
            base_name = os.path.splitext(sFilename_out)[0]
            # Delete all shapefile components
            for ext in ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.qpj', '.sbn', '.sbx']:
                file_to_delete = base_name + ext
                if os.path.exists(file_to_delete):
                    os.remove(file_to_delete)
        else:
            os.remove(sFilename_out)

    pDataset_out = pDriver_out.CreateDataSource(sFilename_out)
    pLayer_out = pDataset_out.CreateLayer('land_ocean_mask', geom_type=ogr.wkbPolygon)

    # Add fields to the output layer
    pField_id = ogr.FieldDefn('id', ogr.OFTInteger)
    pLayer_out.CreateField(pField_id)
    pField_part_type = ogr.FieldDefn('part_type', ogr.OFTString)
    pLayer_out.CreateField(pField_part_type)

    # Start processing geometries
    iCount = 0

    # Set default folder if not provided and create directory only if workspace is specified
    if sWorkspace_out is not None:
        Path(sWorkspace_out).mkdir(parents=True, exist_ok=True)

    for land_geometry in land_geometries:
        pGeometry_mesh = ogr.CreateGeometryFromWkb(land_geometry.wkb)
        envelope = pGeometry_mesh.GetEnvelope()  # Returns a tuple (minX, maxX, minY, maxY)

        # Check whether it is a multipolygon
        sGeometry_name = pGeometry_mesh.GetGeometryName()

        if sGeometry_name == 'POLYGON':
            # Skip Antarctica region
            if envelope[2] < -60:
                continue

            # Add to output layer
            pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
            pFeature_out.SetGeometry(pGeometry_mesh)
            pFeature_out.SetField('id', iCount)
            pFeature_out.SetField('part_type', 'polygon')
            pLayer_out.CreateFeature(pFeature_out)
            pFeature_out.Destroy()

            # Only save individual parts if workspace is specified
            if sWorkspace_out is not None:
                sFilename_part = os.path.join(sWorkspace_out, f'land_geometry_{iCount}.geojson')
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

        elif sGeometry_name == 'MULTIPOLYGON':
            nPart1 = pGeometry_mesh.GetGeometryCount()
            print(nPart1)
            for iPart in range(0, nPart1, 1):
                pGeometry_part = pGeometry_mesh.GetGeometryRef(iPart)
                envelope = pGeometry_part.GetEnvelope()
                if envelope[2] < -60:
                    continue

                # Add to output layer
                pFeature_out = ogr.Feature(pLayer_out.GetLayerDefn())
                pFeature_out.SetGeometry(pGeometry_part)
                pFeature_out.SetField('id', iCount)
                pFeature_out.SetField('part_type', 'multipolygon_part')
                pLayer_out.CreateFeature(pFeature_out)
                pFeature_out.Destroy()

                # Only save individual parts if workspace is specified
                if sWorkspace_out is not None:
                    sFilename_part = os.path.join(sWorkspace_out, f'land_geometry_{iCount}.geojson')
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

    # Clean up output dataset
    pDataset_out.Destroy()

    print(f'Land ocean mask created with {iCount} individual features.')
    return

