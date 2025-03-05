import os
import datetime
import textwrap
import numpy as np
from urllib.error import URLError
from osgeo import  osr, gdal, ogr
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.cm as cm
from matplotlib.collections import PathCollection
from matplotlib.path import Path
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import cartopy as cpl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.toolbox.math.stat.remap import remap
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.map_servers import calculate_zoom_level, calculate_scale_denominator
from pyearth.visual.map.map_servers import StadiaStamen, EsriTerrain, EsriRelief, EsriHydro
from pyearth.visual.map.map_servers import Stadia_terrain_images, Esri_terrain_images, Esri_relief_images, Esri_hydro_images

#get the current year for openstreetmap copy right label

iYear_current = datetime.datetime.now().year
sYear = str(iYear_current)

def map_multiple_vector_files(aFiletype_in,
                             aFilename_in,
                             iFlag_colorbar_in = None,
                             iFlag_title_in = None,
                             iFlag_zebra_in=None,
                             iFlag_filter_in = None,
                             aFlag_thickness_in = None,
                             aFlag_color_in = None,
                             aFlag_discrete_in = None,
                             aFlag_fill_in = None,
                             aVariable_in = None,
                             sFilename_output_in=None,
                             iFlag_scientific_notation_colorbar_in = None,
                             iFlag_openstreetmap_in = None,
                             iFlag_terrain_image_in = None,
                             iFlag_esri_hydro_image_in = None,
                             iFont_size_in=None,
                             iBasemap_zoom_level_in = None,
                             sColormap_in = None,
                             sTitle_in = None,
                             iDPI_in = None,
                             iSize_x_in = None,
                             iSize_y_in = None,
                             aThickness_in = None,
                             aColor_in = None,
                             aMissing_value_in = None,
                             aData_max_in = None,
                             aData_min_in = None,
                             sExtend_in =None,
                             sFont_in = None,
                             sUnit_in=None,
                             aLegend_in = None,
                             aExtent_in = None,
                             pProjection_map_in=None,
                             pProjection_data_in=None,):
    """
    plot vector data on a map
    currently only support geojson and shapefile
    by default, the program will plot all the polygons in the file
    in the furture, the program will support to plot only a subset of polygons

    Because of the overlay effect, it is ideal to plot them in the following order: polygon->polyling->point

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

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    #check vector type first
    if aFiletype_in is None:
        print('Error: please specify the vector type')
        return
    else:
        #point: 1, polyline: 2, polygon: 3
        arr = np.array(aFiletype_in)
        # Check if the array is descending or has the same values
        is_descending = np.all(np.diff(arr) <= 0)
        if is_descending == True:
            pass
        else:
            print('Error: the vector type is not correct')
            return

    pDriver = ogr.GetDriverByName('GeoJSON')
    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSRS_geodetic = ccrs.Geodetic()


    nFile = len(aFilename_in)
    if aFlag_thickness_in is None:
        aFlag_thickness= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_thickness = aFlag_thickness_in

    if aFlag_color_in is None:
        aFlag_color= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_color = aFlag_color_in

    if aFlag_discrete_in is None:
        aFlag_discrete= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_discrete = aFlag_discrete_in

    if aFlag_fill_in is None:
        aFlag_fill= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_fill = aFlag_fill_in

    #get the extent first
    sFilename_in = aFilename_in[0]

    pDataset = ogr.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if iFlag_filter_in is not None:
        iFlag_filter = iFlag_filter_in
    else:
        iFlag_filter = 0

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8

    if iFlag_title_in is not None:
        iFlag_title = iFlag_title_in
        if iFlag_title == 1:
            if sTitle_in is not None:
                sTitle = sTitle_in
            else:
                sTitle =  ''

        else:
            iFlag_title=0

    else:
        iFlag_title = 0
        sTitle =  ''


    if aMissing_value_in is not None:
        aMissing_value = aMissing_value_in
    else:
        aMissing_value =  np.full(nFile, -9999)
    if aData_min_in is not None:
        aData_min = aData_min_in
        aFlag_data_min = np.full(nFile, 1)
    else:
        aData_min=np.full(nFile, -9999)
        aFlag_data_min = np.zeros(nFile, dtype=np.int16)

    if aData_max_in is not None:
        aData_max = aData_max_in
        aFlag_data_max = np.full(nFile, 1)
    else:
        aData_max = np.full(nFile, -9999)
        aFlag_data_max = np.zeros(nFile, dtype=np.int16)
        pass

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'rainbow'

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend =  'max'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit =  ''

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    if aVariable_in is None:
        aVariable = list()
        for i in range(nFile):
            aVariable.append('')
    else:
        aVariable = aVariable_in

    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    cmap = plt.colormaps[sColormap]
    fig = plt.figure( dpi = iDPI  )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    plot_width_inch = fig.get_size_inches()[0] * fig.dpi
    char_width_inch = 0.1 * fig.dpi
    cwidth = int(plot_width_inch / char_width_inch)

    #we require that the first polygon file defines the extent
    pLayer = pDataset.GetLayer(0)
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type =='MULTIPOLYGON':
            for j in range(pGeometry_in.GetGeometryCount()):
                pPolygon = pGeometry_in.GetGeometryRef(j)
                aCoords_gcs = get_geometry_coordinates(pPolygon)
                dLon_max = float(np.max( [dLon_max, np.max(aCoords_gcs[:,0])] ))
                dLon_min = float(np.min( [dLon_min, np.min(aCoords_gcs[:,0])] ))
                dLat_max = float(np.max( [dLat_max, np.max(aCoords_gcs[:,1])] ))
                dLat_min = float(np.min( [dLat_min, np.min(aCoords_gcs[:,1])] ))
        else:
            if sGeometry_type =='POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                dLon_max = float(np.max( [dLon_max, np.max(aCoords_gcs[:,0])] ))
                dLon_min = float(np.min( [dLon_min, np.min(aCoords_gcs[:,0])] ))
                dLat_max = float(np.max( [dLat_max, np.max(aCoords_gcs[:,1])] ))
                dLat_min = float(np.min( [dLat_min, np.min(aCoords_gcs[:,1])] ))
            else:
                if sGeometry_type =='LINESTRING':
                    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    dLon_max = float(np.max( [dLon_max, np.max(aCoords_gcs[:,0])] ))
                    dLon_min = float(np.min( [dLon_min, np.min(aCoords_gcs[:,0])] ))
                    dLat_max = float(np.max( [dLat_max, np.max(aCoords_gcs[:,1])] ))
                    dLat_min = float(np.min( [dLat_min, np.min(aCoords_gcs[:,1])] ))

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                            central_latitude = 0.50*(dLat_max+dLat_min),
                                            globe=None)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection= pProjection_map  ) #projection=ccrs.PlateCarree()
    ax.set_global()
    ax.coastlines(color='black', linewidth=1,resolution='10m')
    # Create an OSM image tile source
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
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    minx, maxx, miny,  maxy = aExtent

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
            print('Basemap zoom level: ',iBasemap_zoom_level)
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
                    esri_terrain = EsriTerrain()
                    esri_hydro = EsriHydro()
                    ll_target_domain = sgeom.box(minx, miny, maxx, maxy)
                    multi_poly = esri_hydro.crs.project_geometry(ll_target_domain, ccrs.PlateCarree())
                    target_domain = multi_poly.geoms[0]
                    _, aExtent_terrain, _ = esri_terrain.image_for_domain(target_domain, iBasemap_zoom_level)
                    _, aExtent_hydro, _ = esri_hydro.image_for_domain(target_domain, iBasemap_zoom_level)
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



    #====================================
    #should we allow more than one scale for one variable?
    aValue_all = list()
    for i in range(nFile):
        aValue = list()
        sFilename = aFilename_in[i]
        iFlag_thickness = aFlag_thickness[i]
        iFlag_color = aFlag_color[i]
        sVariable = aVariable[i]
        pDataset = ogr.Open(sFilename, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        if sVariable is None:
            sVariable=''
        if len(sVariable) > 0:
            #select field value using sqlite
            pKwargs = 'SELECT ' + sVariable + ' FROM ' + pLayer.GetName()
            if iFlag_thickness ==1 :
                pLayer_temp = pDataset.ExecuteSQL(pKwargs)
            else:
                if iFlag_color  == 1:
                    pLayer_temp = pDataset.ExecuteSQL(pKwargs)
            values = []
            # Iterate over the layer
            for feature in pLayer_temp:
                # Get the desired field value
                value = feature.GetField(sVariable)
                # Append the value to the list
                values.append(value)
            #conver aValue to numpy array
            aValue = np.array(values)
            aValue_all.append(aValue)
        else:
            aValue_all.append(aValue)

    iThickness_max = 2.5
    iThickness_min = 0.3

    #set min and max
    for i in range(nFile):
        iFlag_data_min = aFlag_data_min[i]
        iFlag_data_max = aFlag_data_max[i]
        dMissing_value = aMissing_value[i]
        aValue = np.array(aValue_all[i])
        aValue = aValue[aValue != dMissing_value]
        if len(aValue) == 0:
            aData_min[i] = 0
            aData_max[i] = 0
        else:
            if iFlag_data_min != 1:
                aData_min[i] = np.min(aValue)
            if iFlag_data_max != 1:
                aData_max[i] = np.max(aValue)

    for i in range(nFile):
        sFilename = aFilename_in[i]
        iFlag_thickness = aFlag_thickness[i]
        iFlag_color = aFlag_color[i]
        iFlag_discrete = aFlag_discrete[i]
        iFlag_fill = aFlag_fill[i]
        dValue_min = aData_min[i]
        dValue_max = aData_max[i]
        dMissing_value = aMissing_value[i]
        sVariable = aVariable[i]
        aValue = np.array(aValue_all[i])
        pDataset = ogr.Open(sFilename, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        nValue = pLayer.GetFeatureCount()
        if iFlag_discrete == 1: #discrete
            aValue = np.array(aValue, dtype=int)
            #reorder it
            aValue = np.unique(aValue)
            aValue = np.sort(aValue)
            nValue_discrete = len(aValue)
            aIndex = np.linspace(0,1,nValue_discrete)
            prng = np.random.RandomState(1234567890)
            prng.shuffle(aIndex)
            colors = plt.colormaps[sColormap](aIndex)
            aColorMap = ListedColormap(colors)
            pass
        else: #continuous
            aColorMap = plt.colormaps[sColormap](np.linspace(0, 1, nValue))
        lID = 0
        #aPoint=list()
        aPoint_x=list()
        aPoint_y=list()
        aMarker = list()
        aSize=list()
        aPolyline=list()
        aThickness=list()
        aPolygon = list()
        aColor = list()

        if iFlag_filter == 1:
            pLayer.SetSpatialFilterRect(minx, maxx, miny, maxy)

        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if aThickness_in is not None:
                iThickness = aThickness_in[i]
            else:
                if iFlag_thickness ==1 :
                    dValue = float(pFeature.GetField(sVariable))
                    if dValue < dValue_min or dValue > dValue_max:
                        continue

                    dValue = np.clip(dValue, dValue_min, dValue_max)

                    iThickness = remap( dValue, dValue_min, dValue_max, iThickness_min, iThickness_max )
                else:
                    iThickness = 0.25

            if aColor_in is not None:
                sColor = aColor_in[i]
            else:
                if iFlag_color == 1:
                    dValue = float(pFeature.GetField(sVariable))
                    if iFlag_discrete ==1:
                        if dValue < dValue_min or dValue > dValue_max:
                            continue
                        iValue = int(dValue)
                        iColor_index = np.where(aValue == iValue)[0][0]
                        sColor = aColorMap(iColor_index)
                    else:
                        if dValue < dValue_min or dValue > dValue_max:
                            continue

                        iColor_index = int( (dValue - dValue_min) / (dValue_max - dValue_min) * 255 )
                        sColor = aColorMap[iColor_index]
                else:
                    sColor = 'black'

            if sGeometry_type =='POINT':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = aCoords_gcs[:,0:2]
                aColor.append(sColor)
                aPoint_x.append(aCoords_gcs[0,0])
                aPoint_y.append(aCoords_gcs[0,1])
                aSize.append(iThickness)
                aMarker.append('o')
            else:
                if sGeometry_type =='LINESTRING':
                    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    aCoords_gcs = aCoords_gcs[:,0:2]
                    nvertex = len(aCoords_gcs)
                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpl.path.Path.MOVETO
                    path = mpl.path.Path(aCoords_gcs, codes)
                    x, y = zip(*path.vertices)
                    aThickness.append(iThickness)
                    aColor.append(sColor)
                    aPolyline.append(list(zip(x, y)))
                else:
                    if sGeometry_type == 'POLYGON':
                        aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                        aColor.append(sColor)
                        aPolygon.append(aCoords_gcs[:, 0:2])
                    else:
                        if sGeometry_type == 'MULTIPOLYGON':
                            for j in range(pGeometry_in.GetGeometryCount()):
                                pPolygon = pGeometry_in.GetGeometryRef(j)
                                aCoords_gcs = get_geometry_coordinates(pPolygon)
                                aColor.append(sColor)
                                aPolygon.append(aCoords_gcs[:, 0:2])


            lID = lID + 1

        if len(aPolygon) > 0:
            aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
            if iFlag_fill == 1:
                pPC = PatchCollection(aPatch, cmap=cmap, alpha=0.8, edgecolor=None,
                                      facecolor=aColor, linewidths=0.25,
                                      transform=pProjection_data)
            else:
                pPC = PatchCollection(aPatch, cmap=cmap, alpha=0.8, edgecolor=aColor,
                                      facecolor='none', linewidths=0.25,
                                      transform=pProjection_data)
            ax.add_collection(pPC)

        if len(aPolyline) > 0:
            #polyline
            pLC = LineCollection(aPolyline,  alpha=0.8, edgecolor=aColor,
                         facecolor='none', linewidths=aThickness, transform=pProjection_data)
            ax.add_collection(pLC)

        if len(aPoint_x) > 0:
            for i in range(len(aPoint_x)):
                ax.scatter(aPoint_x[i], aPoint_y[i], c=aColor[i],
                                 s=aSize[i], marker=aMarker[i],
                                   alpha=dAlpha, transform=pProjection_data)

    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
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
        ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
        ax.set_axis_off()

    #reset extent
    ax.set_extent(aExtent, crs = pSRS_wgs84)

    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass
        ax.set_title(sTitle)

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            #dLocation = 0.06 + i * 0.04
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size)


    if iFlag_colorbar ==1:
        #ax_cb= fig.add_axes([0.75, 0.15, 0.02, 0.6])
        fig.canvas.draw()
        ax_pos = ax.get_position() # get the original position
        ax_cb = fig.add_axes([ax_pos.x1+0.06, ax_pos.y0, 0.02, ax_pos.height])
        if iFlag_scientific_notation_colorbar==1:
            formatter = OOMFormatter(fformat= "%1.1e")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=cmap,
                                           norm=mpl.colors.Normalize(dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
        else:
            formatter = OOMFormatter(fformat= "%1.1f")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=cmap,
                                           norm=mpl.colors.Normalize(dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 10
        cb.ax.set_ylabel(sUnit, rotation=270)
        cb.ax.tick_params(labelsize=6)

    if iFlag_zebra ==1:
        ax.set_axis_off()
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent, crs = pSRS_wgs84)
    pDataset = pLayer = pFeature  = None

    if sFilename_output_in is None:
        plt.show()
    else:
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
                plt.savefig(sFilename_out, bbox_inches='tight', format ='ps')

        #clean cache
        plt.close('all')
        plt.clf()
