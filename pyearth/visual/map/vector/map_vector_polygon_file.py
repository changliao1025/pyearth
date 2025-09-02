import os
import datetime
import textwrap
import numpy as np
from osgeo import osr, gdal, ogr
from urllib.error import URLError
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as mpolygon
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyproj import Geod
import shapely.geometry as sgeom
from shapely.geometry import Polygon as spolygon
from pyearth.system.define_global_variables import *
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.gis.geometry.calculate_polygon_area import calculate_polygon_area
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.map.map_servers import StadiaStamen, EsriTerrain, EsriRelief, EsriHydro
from pyearth.visual.map.map_servers import Stadia_terrain_images, Esri_terrain_images, Esri_relief_images, Esri_hydro_images

osr.UseExceptions()

#get the current year for openstreetmap copy right label
iYear_current = datetime.datetime.now().year
sYear = str(iYear_current)
def simplify_coordinates(aCoords_gcs, min_distance=1.0):
    """
    Simplify a set of geographic coordinates by retaining points that are at least
    a specified minimum distance apart.

    Args:
        aCoords_gcs (numpy.ndarray): Array of geographic coordinates (longitude, latitude).
        min_distance (float): Minimum distance (in the same units as the input coordinates)
                              between consecutive points to retain.

    Returns:
        numpy.ndarray: Simplified array of geographic coordinates.
    """
    # Initialize list for simplified coordinates
    aCoords_gcs_simplified = []

    # Iterate through points and calculate distances
    nPoint = len(aCoords_gcs)
    for i in range(nPoint - 1):
        dLon, dLat = aCoords_gcs[i]
        dLon1, dLat1 = aCoords_gcs[i + 1]
        dDistance = calculate_distance_based_on_longitude_latitude(dLon, dLat, dLon1, dLat1)
        if dDistance >= min_distance:
            aCoords_gcs_simplified.append((dLon1, dLat1))

    # Close the ring by appending the first point
    if aCoords_gcs_simplified:
        aCoords_gcs_simplified.append(aCoords_gcs_simplified[0])

    # Convert to numpy array
    return np.array(aCoords_gcs_simplified)

def map_vector_polygon_file(iFiletype_in,
                            sFilename_in,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFlag_color_in = None,
                            iFlag_colorbar_in=None,
                            iFlag_zebra_in=None,
                            iFlag_fill_in=None,
                            iFont_size_in=None,
                            iFlag_discrete_in=None,
                            iFlag_filter_in = None,
                            iFlag_openstreetmap_in = None,
                            iFlag_terrain_image_in = None,
                            iFlag_relief_image_in = None,
                            iFlag_esri_hydro_image_in = None,
                            iFlag_buffer_in = None,
                            iFlag_debug_in = None,
                            iBasemap_zoom_level_in = None,
                            sColor_in = None,
                            sColormap_in=None,
                            sTitle_in=None,
                            iDPI_in=None,
                            iSize_x_in=None,
                            iSize_y_in=None,
                            dMissing_value_in=None,
                            dData_max_in=None,
                            dData_min_in=None,
                            dLinewidth_in=None,
                            sField_color_in = None,
                            sExtend_in=None,
                            sFont_in=None,
                            sUnit_in=None,
                            sFormat_colorbar_in=None,
                            aLegend_in=None,
                            aExtent_in=None,
                            aDict_discrete_labels_in = None,
                            aDict_value_color_in = None,
                            pProjection_map_in=None,
                            pProjection_data_in = None):
    """
    plot vector data on a map
    currently only support geojson and shapefile
    by default, the program will plot all the polygons in the file
    in the furture, the program will support to plot only a subset of polygons

    Args:
        iFiletype_in (_type_): _description_
        sFilename_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
    """
    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSRS_geodetic = ccrs.Geodetic()

    if iFiletype_in == 1:  # geojson
        pDriver = ogr.GetDriverByName('GeoJSON')
    else:
        if iFiletype_in == 2:  # shapefile
            pDriver = ogr.GetDriverByName('Esri Shapefile')
        else:
            pDriver = ogr.GetDriverByName('Parquet')

    if pDriver is None:
        print('Driver not available')
        return

    if os.path.exists(sFilename_in) is False:
        print('File does not exist')
        return

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    if pDataset is None:
        print('Dataset not available')
        return

    pLayer = pDataset.GetLayer(0)
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 150

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if iFlag_colorbar ==1:
        iFlag_color = 1

    if iFlag_fill_in is not None:
        iFlag_fill = iFlag_fill_in
    else:
        iFlag_fill = True

    if iFlag_discrete_in is not None:
        iFlag_discrete = iFlag_discrete_in
    else:
        iFlag_discrete = 0

    if iFlag_filter_in is not None:
        iFlag_filter = iFlag_filter_in
    else:
        iFlag_filter = 0

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

    if dData_min_in is not None:
        iFlag_data_min = 1
        dData_min = dData_min_in
    else:
        iFlag_data_min = 0
        pass

    if dData_max_in is not None:
        iFlag_data_max = 1
        dData_max = dData_max_in
    else:
        iFlag_data_max = 0
        pass

    if dLinewidth_in is not None:
        dLinewidth = dLinewidth_in
    else:
        dLinewidth = 0.25

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap = 'rainbow'

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title = 1
    else:
        iFlag_title = 0
        sTitle = ''
    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend = 'max'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit = ''

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    if sFilename_output_in is None:
        plt.ion()

    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aValue = list()

    if iFlag_color == 1:
        #check color field
        if sField_color_in  is not None:
            sField_color = sField_color_in
        else:
            print('Color field is not provided')
            return

    nFeature = pLayer.GetFeatureCount()

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_color == 1:
            dValue = float(pFeature.GetField(sField_color))
            aValue.append(dValue)

        if sGeometry_type == 'MULTIPOLYGON':
            for i in range(pGeometry_in.GetGeometryCount()):
                pPolygon = pGeometry_in.GetGeometryRef(i)
                aCoords_gcs = get_geometry_coordinates(pPolygon)
                dLon_max = float(np.max([dLon_max, np.max(aCoords_gcs[:, 0])]))
                dLon_min = float(np.min([dLon_min, np.min(aCoords_gcs[:, 0])]))
                dLat_max = float(np.max([dLat_max, np.max(aCoords_gcs[:, 1])]))
                dLat_min = float(np.min([dLat_min, np.min(aCoords_gcs[:, 1])]))
        else:
            if sGeometry_type == 'POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                dLon_max = float(np.max([dLon_max, np.max(aCoords_gcs[:, 0])]))
                dLon_min = float(np.min([dLon_min, np.min(aCoords_gcs[:, 0])]))
                dLat_max = float(np.max([dLat_max, np.max(aCoords_gcs[:, 1])]))
                dLat_min = float(np.min([dLat_min, np.min(aCoords_gcs[:, 1])]))

    if iFlag_color == 1:
        aValue = np.array(aValue)
        if iFlag_discrete ==1:
            #convert to integer
            aValue = np.array(aValue, dtype=int)
            aValue = np.unique(aValue)
            aValue = np.sort(aValue)
            nValue = len(aValue) #get unique values from the aValue

        if iFlag_data_min == 1:  # min is provided
            dValue_min = dData_min  # np.min(aValue)
        else:
            aValue = aValue[aValue != dMissing_value]
            dValue_min = np.min(aValue)

        if iFlag_data_max == 1:  # max is provided
            dValue_max = dData_max  # np.max(aValue)
        else:
            aValue = aValue[aValue != dMissing_value]
            dValue_max = np.max(aValue)

        aValue = np.clip(aValue, a_min=dValue_min, a_max=dValue_max)

        print(dValue_min, dValue_max)
        if dValue_max == dValue_min:
            return

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        dLon_mean = 0.50 * (dLon_max + dLon_min)
        dLat_mean = 0.50 * (dLat_max + dLat_min)
        pProjection_map = ccrs.Orthographic(central_longitude=dLon_mean, central_latitude=dLat_mean)
        print(dLon_mean, dLat_mean)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84

    pProjection_map._threshold /= 1.0E6

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map)
    plot_width_inch = fig.get_size_inches()[0] * fig.dpi
    char_width_inch = 0.1 * fig.dpi
    cwidth = int(plot_width_inch / char_width_inch)

    if iFlag_discrete == 1:
        nValue_current = len(aValue)
        if aDict_value_color_in is not None:
            custom_colors = []
            all_values_have_colors = True
            for val_code in aValue: # aValue is sorted unique values
                color_for_val = aDict_value_color_in.get(val_code)
                if color_for_val is None:
                    print(f"Warning: Color for value {val_code} not found in dict_value_colors_in. Falling back to sColormap.")
                    all_values_have_colors = False
                    break
                custom_colors.append(color_for_val)
            if all_values_have_colors:
                pCmap = ListedColormap(custom_colors)
            else: # Fallback if custom color mapping is incomplete
                if nValue_current > 0:
                    # Fallback to using sColormap (consider if shuffling is always desired)
                    aIndex = np.linspace(0, 1, nValue_current)
                    # prng = np.random.RandomState(1234567890) # Shuffling might not be ideal for predefined qualitative maps
                    # prng.shuffle(aIndex)
                    colors_from_cmap = plt.colormaps[sColormap](aIndex)
                    pCmap = ListedColormap(colors_from_cmap)
                else:
                    pCmap = ListedColormap(['white']) # Default for no values
        else:
            aIndex = np.linspace(0,1,nValue)
            prng = np.random.RandomState(1234567890)
            prng.shuffle(aIndex)
            #print(aIndex)
            colors = plt.colormaps[sColormap](aIndex)
            pCmap = ListedColormap(colors)
    else:
        pCmap = plt.colormaps[sColormap]

    if aExtent_in is None:
        marginx = (dLon_max - dLon_min) / 20
        marginy = (dLat_max - dLat_min) / 20
        if (dLat_max + marginy)> 90:
            dLat_max = 90
        else:
            dLat_max = dLat_max + marginy
        if (dLat_min - marginy) < -90:
            dLat_min = -90
        else:
            dLat_min = dLat_min - marginy
        if (dLon_max + marginx) > 180:
            dLon_max = 180
        else:
            dLon_max = dLon_max + marginx
        if (dLon_min - marginx) < -180:
            dLon_min = -180
        else:
            dLon_min = dLon_min - marginx
        aExtent = [dLon_min, dLon_max, dLat_min, dLat_max]
    else:
        aExtent = aExtent_in

    print(aExtent)
    minx, maxx, miny, maxy = aExtent
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    if iFlag_filter == 1:
        pLayer.SetSpatialFilterRect(minx, maxx, miny, maxy)
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    ax.coastlines(linewidth=0.5, color='k', resolution='10m')
    try:
        dAlpha = 1.0
        #only one of the base map can be used
        if iBasemap_zoom_level_in is not None:
            iBasemap_zoom_level = iBasemap_zoom_level_in
        else:
            image_size = [1000, 1000]
            scale_denominator = calculate_scale_denominator(aExtent, image_size)
            pSrc = osr.SpatialReference()
            pSrc.ImportFromEPSG(3857) # mercator
            pProjection = pSrc.ExportToWkt()
            iBasemap_zoom_level = calculate_zoom_level(scale_denominator, pProjection, dpi=int(iDPI/2))
            #print('Basemap zoom level: ',iBasemap_zoom_level)
            pass
        if iFlag_openstreetmap_in is not None and iFlag_openstreetmap_in == 1:
            from cartopy.io.img_tiles import OSM
            osm_tiles = OSM()
            #Add the OSM image to the map
            ax.add_image(osm_tiles, iBasemap_zoom_level) #, alpha=0.5
            sLicense_info = "© OpenStreetMap contributors "+ sYear + "." + " Distributed under the Open Data Commons Open Database License (ODbL) v1.0."

            sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=cwidth))
            ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                    color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

            #we also need to set transparency for the image to be added
            dAlpha = 0.5
        else:
            if iFlag_terrain_image_in is not None and iFlag_terrain_image_in == 1:
                stamen_terrain = StadiaStamen()
                ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
                multi_poly = stamen_terrain.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
                target_domain = multi_poly.geoms[0]
                _, aExtent_terrain, _ = stamen_terrain.image_for_domain(target_domain, iBasemap_zoom_level)
                img_stadia_terrain = Stadia_terrain_images(aExtent, iBasemap_zoom_level)
                ax.imshow(img_stadia_terrain,  extent=aExtent_terrain, transform=stamen_terrain.crs)
                #add the license information
                sLicense_info = "© Stamen Design, under a Creative Commons Attribution (CC BY 3.0) license."
                sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=cwidth))
                ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                        color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
                dAlpha = 0.8
            else:
                if iFlag_esri_hydro_image_in is not None and iFlag_esri_hydro_image_in == 1:
                    if iFlag_relief_image_in is not None:
                        esri_terrain = EsriRelief()
                    else:
                        esri_terrain = EsriTerrain()
                    esri_hydro = EsriHydro()
                    ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
                    multi_poly = esri_hydro.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
                    target_domain = multi_poly.geoms[0]
                    _, aExtent_terrain, _ = esri_terrain.image_for_domain(target_domain, iBasemap_zoom_level)
                    _, aExtent_hydro, _ = esri_hydro.image_for_domain(target_domain, iBasemap_zoom_level)
                    if iFlag_relief_image_in is not None:
                        img_esri_terrain = Esri_relief_images(aExtent, iBasemap_zoom_level)
                    else:
                        img_esri_terrain  = Esri_terrain_images(aExtent, iBasemap_zoom_level)
                    img_eari_hydro  = Esri_hydro_images(aExtent, iBasemap_zoom_level)
                    ax.imshow(img_esri_terrain,  extent=aExtent_terrain, transform=esri_terrain.crs)
                    ax.imshow(img_eari_hydro,  extent=aExtent_hydro, transform=esri_hydro.crs)
                    #add the license information
                    sLicense_info = "© Esri Hydro Reference Overlay"
                    sLicense_info_wrapped = "\n".join(textwrap.wrap(sLicense_info, width=60))
                    ax.text(0.5, 0.05, sLicense_info_wrapped, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                            color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
                    dAlpha = 0.8
                else:
                    pass

    except URLError as e:
        print('No internet connection')
        dAlpha = 1.0

    aPolygon = list()
    aColor_index=list()
    aEdgecolor = list()
    #aFacecolor = list()
    if iFlag_debug_in is not None:
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            # get attribute
            if iFlag_color == 1:
                dValue = float(pFeature.GetField(sField_color))
                if dValue != dMissing_value:
                    if dValue > dValue_max:
                        dValue = dValue_max
                        #continue
                    if dValue < dValue_min:
                        dValue = dValue_min
                        #continue
                    if iFlag_discrete ==1:
                        #use unique value method to assign color
                        iValue = int(dValue) #this maight cause issue if the date is not using integer
                        iColor_index = np.where(aValue == iValue)[0][0]
                    else:
                        iColor_index = int((dValue - dValue_min) /
                                       (dValue_max - dValue_min) * 255)
                else:
                    iColor_index = 0
                    pass
            else:
                iColor_index = 0
                pass
            if sGeometry_type == 'MULTIPOLYGON':
                for i in range(pGeometry_in.GetGeometryCount()):
                    pPolygon = pGeometry_in.GetGeometryRef(i)
                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                    aCoords_gcs=aCoords_gcs[:,0:2]
                    aColor_index.append(iColor_index)
                    aEdgecolor.append('none')
                    #manually plot the polygon
                    mPolygon = mpolygon(aCoords_gcs, closed=True, edgecolor='blue', facecolor='red',
                                      alpha=dAlpha,  transform=pProjection_data)
                    # Add the polygon to the map
                    ax.add_patch(mPolygon)
            else:
                if sGeometry_type == 'POLYGON':
                    if iFlag_buffer_in is not None:
                        aCoords_gcs = get_geometry_coordinates(pGeometry_in.Buffer(1.0E-5))
                    else:
                        aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    aCoords_gcs = aCoords_gcs[:,0:2]
                    dArea = calculate_polygon_area(aCoords_gcs[:,0], aCoords_gcs[:,1])
                    if dArea > 1.0:
                        if iFlag_fill == 1:
                            mPolygon = mpolygon(aCoords_gcs, closed=True, edgecolor='blue', facecolor='red',
                                      alpha=dAlpha, transform=pProjection_data)
                        else:
                            mPolygon = mpolygon(aCoords_gcs, closed=True, edgecolor='blue', facecolor='none',
                                      alpha=dAlpha, transform=pProjection_data)
                        ax.add_patch(mPolygon)
                    else:
                        pass
                else:
                    print('Geometry type not supported')
    else:
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            # get attribute
            if iFlag_color == 1:
                dValue = float(pFeature.GetField(sField_color))
                if dValue != dMissing_value:
                    if dValue > dValue_max:
                        dValue = dValue_max
                        #continue
                    if dValue < dValue_min:
                        dValue = dValue_min
                        #continue
                    if iFlag_discrete ==1:
                        #use unique value method to assign color
                        iValue = int(dValue) #this maight cause issue if the date is not using integer
                        iColor_index = np.where(aValue == iValue)[0][0]
                    else:
                        iColor_index = int((dValue - dValue_min) /
                                       (dValue_max - dValue_min) * 255)
                else:
                    iColor_index = 0
                    pass
            else:
                iColor_index = 0
                pass
            if sGeometry_type == 'MULTIPOLYGON':
                for i in range(pGeometry_in.GetGeometryCount()):
                    pPolygon = pGeometry_in.GetGeometryRef(i)
                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                    aCoords_gcs=aCoords_gcs[:,0:2]
                    aColor_index.append(iColor_index)
                    aPolygon.append(aCoords_gcs)
                    aEdgecolor.append('none')

            else:
                if sGeometry_type == 'POLYGON':
                    if iFlag_buffer_in is not None:
                        aCoords_gcs = get_geometry_coordinates(pGeometry_in.Buffer(1.0E-5))
                    else:
                        aCoords_gcs = get_geometry_coordinates(pGeometry_in)

                    aCoords_gcs=aCoords_gcs[:,0:2]
                    aColor_index.append(iColor_index)
                    aPolygon.append(aCoords_gcs)
                    aEdgecolor.append('none')
                else:
                    print('Geometry type not supported')

        aColor_index = np.array(aColor_index)
        #flatten the array as 1D
        aColor_index= aColor_index.flatten()
        aColor= pCmap(aColor_index)
        if iFlag_color == 1:
            if iFlag_fill == True:
                aPatch = [mpolygon(poly, closed=True) for poly in aPolygon]
                pPC = PatchCollection(aPatch, alpha=dAlpha,
                                      edgecolor=aColor,
                                      facecolor=aColor,
                                      transform=pProjection_data)
            else:
                aPatch = [mpolygon(poly, closed=True) for poly in aPolygon]
                pPC = PatchCollection(aPatch, alpha=dAlpha,
                                      edgecolor=aColor,
                                      facecolor='none',
                                      linewidths=dLinewidth,
                                      transform=pProjection_data)
        else:
            if sColor_in is not None:
                sColor = sColor_in
            else:
                sColor = 'blue'
            if iFlag_fill == True:
                aPatch = [mpolygon(poly, closed=True, fill=iFlag_fill) for poly in aPolygon]
                pPC = PatchCollection(aPatch, alpha=dAlpha,
                                      edgecolor='none',
                                      facecolor=sColor,
                                      transform=pProjection_data)
            else:
                aPatch = [mpolygon(poly, closed=True, fill=iFlag_fill) for poly in aPolygon]
                pPC = PatchCollection(aPatch, alpha=dAlpha,
                                      edgecolor=sColor,
                                      facecolor='none',
                                      linewidths=dLinewidth,
                                      transform=pProjection_data)

        ax.add_collection(pPC)

    ax.set_extent(aExtent, crs = pSRS_wgs84)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    if iFlag_zebra == 1:
        ax.set_xticks(np.arange(minx, maxx+(maxx-minx)/11, (maxx-minx)/10))
        dummy = (maxy-miny)/10
        ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
        ax.set_axis_off()

    sTitle = "\n".join(textwrap.wrap(sTitle, width=cwidth))
    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.05, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size + 2)


    if iFlag_zebra ==1:
        #ax.set_axis_off()
        ax.set_extent(aExtent, crs = pSRS_wgs84)
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent, crs = pSRS_wgs84)

    if iFlag_colorbar == 1:
        fig.canvas.draw()
        # Section 2
        ax_pos = ax.get_position() # get the original position
        ax_cb = fig.add_axes([ax_pos.x1+0.06, ax_pos.y0, 0.02, ax_pos.height])
        if iFlag_scientific_notation_colorbar == 1:
            formatter = OOMFormatter(fformat="%1.1e")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
        else:
            if iFlag_discrete ==1:
                if len(aValue) > 1:
                    boundaries = [aValue[0] - 0.5] # Start boundary
                    for k_val in range(len(aValue) - 1):
                        boundaries.append((aValue[k_val] + aValue[k_val+1]) / 2.0)
                    boundaries.append(aValue[-1] + 0.5) # End boundary
                elif len(aValue) == 1: # Handle case with only one unique value
                    boundaries = [aValue[0] - 0.5, aValue[0] + 0.5]
                else: # No values, empty colorbar
                    boundaries = [0, 1]

                norm = mpl.colors.BoundaryNorm(boundaries, pCmap.N, clip=True)

                # For discrete colorbars, 'extend' might not be desired unless explicitly handled
                current_extend = sExtend if not iFlag_discrete else 'neither' # Or None, or handle based on actual data clipping
                # Calculate tick positions to be at the center of each color patch
                tick_centers = []
                if len(aValue) > 0: # Ensure aValue is not empty
                    for k_idx in range(len(aValue)): # Iterate for each category/color
                        tick_centers.append((boundaries[k_idx] + boundaries[k_idx+1]) / 2.0)

                #formatter = mpl.ticker.FuncFormatter(lambda x, pos: "{:.0f}".format(x))
                #cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                #                           cmap=pCmap,
                #                           norm=mpl.colors.Normalize(
                #                               dValue_min, dValue_max),  # vmax and vmin
                #                           extend=sExtend, format=formatter)
                cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                               cmap=pCmap,
                                               norm=norm, # Use the BoundaryNorm
                                               extend=current_extend,
                                               ticks=tick_centers)

                discrete_tick_labelsize = iFont_size - 4
                if aDict_discrete_labels_in is not None and len(aValue) > 0:
                    labels = [aDict_discrete_labels_in.get(int(tick_val), str(int(tick_val))) for tick_val in aValue]
                    cb.set_ticklabels(labels)
                    cb.ax.tick_params(axis='y', which='major', labelsize=discrete_tick_labelsize) # Apply specific size
                elif len(aValue) > 0:
                    cb.set_ticklabels([str(int(val)) for val in aValue])
                    cb.ax.tick_params(axis='y', which='major', labelsize=discrete_tick_labelsize)
            else:
                if sFormat_colorbar_in is not None:
                    formatter = sFormat_colorbar_in
                else:
                    formatter = OOMFormatter(fformat="%1.2f")
                cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 3
        cb.ax.set_ylabel(sUnit, rotation=90, fontsize=iFont_size-2)
        cb.ax.get_yaxis().set_label_position('left')
        cb.ax.tick_params(labelsize=iFont_size-2)


    pDataset = pLayer = pFeature = None


    if sFilename_output_in is None:
        plt.show(block=True)
        print('Finished plotting interactive map')
    else:
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == '.png':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            if sExtension == '.pdf':
                plt.savefig(sFilename_out, bbox_inches='tight')
            else:
                plt.savefig(sFilename_out, bbox_inches='tight', format='ps')

        plt.close('all')
        plt.clf()
        print('The plot is saved to: ', sFilename_out)
