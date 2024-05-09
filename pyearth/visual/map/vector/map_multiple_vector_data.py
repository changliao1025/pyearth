import os
import numpy as np
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
from cartopy.io.img_tiles import OSM
from pyearth.toolbox.math.stat.remap import remap
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter

#get the current year for openstreetmap copy right label
import datetime
iYear_current = datetime.datetime.now().year
#convert to string
sYear = str(iYear_current)

def map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             iFlag_colorbar_in = None,
                             iFlag_title_in = None,
                             aFlag_thickness_in = None,
                             aFlag_color_in = None,
                             aFlag_discrete_in = None,
                             aFlag_fill_in = None,
                             aVariable_in = None,
                             sFilename_output_in=None,
                             iFlag_scientific_notation_colorbar_in = None,
                             iFlag_openstreetmap_in = None,
                             iFont_size_in=None,
                             iFlag_openstreetmap_level_in = None,
                             sColormap_in = None,
                             sTitle_in = None,
                             iDPI_in = None,
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
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

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

    plt.rcParams["font.family"] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    cmap = cm.get_cmap(sColormap)
    fig = plt.figure( dpi = iDPI  )
    iSize_x= 8
    iSize_y= 8
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    #we require that the first polygon file defines the extent
    pLayer = pDataset.GetLayer(0)
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
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
                dLon_max = np.max( [dLon_max, np.max(aCoords_gcs[:,0])] )
                dLon_min = np.min( [dLon_min, np.min(aCoords_gcs[:,0])] )
                dLat_max = np.max( [dLat_max, np.max(aCoords_gcs[:,1])] )
                dLat_min = np.min( [dLat_min, np.min(aCoords_gcs[:,1])] )
        else:
            if sGeometry_type =='POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                dLon_max = np.max( [dLon_max, np.max(aCoords_gcs[:,0])] )
                dLon_min = np.min( [dLon_min, np.min(aCoords_gcs[:,0])] )
                dLat_max = np.max( [dLat_max, np.max(aCoords_gcs[:,1])] )
                dLat_min = np.min( [dLat_min, np.min(aCoords_gcs[:,1])] )
            else:
                if sGeometry_type =='LINESTRING':
                    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    dLon_max = np.max( [dLon_max, np.max(aCoords_gcs[:,0])] )
                    dLon_min = np.min( [dLon_min, np.min(aCoords_gcs[:,0])] )
                    dLat_max = np.max( [dLat_max, np.max(aCoords_gcs[:,1])] )
                    dLat_min = np.min( [dLat_min, np.min(aCoords_gcs[:,1])] )

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

    # Create an OSM image tile source
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in

    ax.set_global()
    print(aExtent)
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    if iFlag_openstreetmap_in is not None and iFlag_openstreetmap_in == 1:
        if iFlag_openstreetmap_level_in is not None:
            iFlag_openstreetmap_level = iFlag_openstreetmap_level_in
        else:
            iFlag_openstreetmap_level = 9
            pass

        osm_tiles = OSM()
        #Add the OSM image to the map
        ax.add_image(osm_tiles, iFlag_openstreetmap_level)
        sLicense_info = "© OpenStreetMap contributors "+ sYear + "." + " Distributed under the Open Data Commons Open Database License (ODbL) v1.0."
        ax.text(0.5, 0.05, sLicense_info, transform=ax.transAxes, ha='center', va='center', fontsize=6,
                color='gray', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))


    minx, miny, maxx, maxy = aExtent
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
        sVariable = aVariable_in[i]
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
            colors = plt.cm.get_cmap(sColormap)(aIndex)
            aColorMap = ListedColormap(colors)
            pass
        else: #continuous
            aColorMap = plt.cm.get_cmap(sColormap)(np.linspace(0, 1, nValue))
        lID = 0
        aPoint=list()
        aPolyline=list()
        aThickness=list()
        aPolygon = list()
        aColor = list()
        pLayer.SetSpatialFilterRect(minx, miny, maxx, maxy)

        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if iFlag_thickness ==1 :
                dValue = float(pFeature.GetField(sVariable))
                if dValue < dValue_min or dValue > dValue_max:
                    continue

                dValue = np.clip(dValue, dValue_min, dValue_max)

                iThickness = remap( dValue, dValue_min, dValue_max, iThickness_min, iThickness_max )
            else:
                iThickness = 0.25

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
                #pick color from colormap
            if sGeometry_type =='POINT':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = aCoords_gcs[:,0:2]
                #ax.plot(aCoords_gcs[0], aCoords_gcs[1], 'o', color= sColor, markersize=2, transform=ccrs.Geodetic())
                aColor.append(sColor)
                aPoint.append(aCoords_gcs)
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

        if len(aPoint) > 0:
            paths = [Path([point]) for point in aPoint]
            pPC = PathCollection(paths, alpha=0.8, edgecolor=aColor,
                                 facecolor=aColor, linewidths=aThickness, transform=pProjection_data)
            ax.add_collection(pPC)

        if len(aPolyline) > 0:
            #polyline
            pLC = LineCollection(aPolyline,  alpha=0.8, edgecolor=aColor,
                         facecolor='none', linewidths=aThickness, transform=pProjection_data)
            ax.add_collection(pLC)

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


    #reset extent
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    ax.coastlines(color='black', linewidth=0.5)
    iFlag_label = 0
    if iFlag_label == 1:
        sText = 'Manaus'
        dLongitude_label = -60.016667
        dLatitude_label  = -3.1
        ax.text(dLongitude_label, dLatitude_label, sText,
                verticalalignment='center', horizontalalignment='center',
                color='black', fontsize=iFont_size,transform=pProjection_data)

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

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k', 'rotation':90,'weight': 'normal'}
    if iFlag_title==1:
        ax.set_title( sTitle )

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
