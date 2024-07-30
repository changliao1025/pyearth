import os
import numpy as np
from osgeo import  osr, gdal, ogr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap

import cartopy as cpl
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter
from pyearth.toolbox.math.stat.remap import remap
from pyearth.visual.map.zebra_frame import zebra_frame

def map_vector_polyline_file(iFiletype_in,
                             sFilename_in,
                             sVariable_in=None,
                             sFilename_output_in=None,
                             iFlag_thickness_in =None,
                             iFlag_color_in= None,
                             iFlag_colorbar_in=None,
                             iFlag_zebra_in=None,
                             iFlag_scientific_notation_colorbar_in=None,
                             iFlag_label_in = None,
                             iFlag_filter_in = None,
                             iFlag_discrete_in=None,
                             iDPI_in = None,
                             iSize_x_in = None,
                             iSize_y_in = None,
                             iFont_size_in = None,
                             sField_thickness_in=None,
                             sField_color_in=None,
                             sColormap_in = None,
                             sTitle_in = None,
                             dMissing_value_in=None,
                             dData_max_in = None,
                             dData_min_in = None,
                             sFont_in = None,
                             sExtend_in =None,
                             sUnit_in=None,
                             aLegend_in = None,
                             aExtent_in = None,
                             pProjection_map_in = None,
                             pProjection_data_in = None):
    """
    Map a vector polyline data

    Args:
        iFiletype_in (int): _description_
        sFilename_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        iFlag_thickness_in (_type_, optional): _description_. Defaults to None.
        sField_thickness_in (_type_, optional): _description_. Defaults to None.
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        sFont_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
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

    if os.path.exists(sFilename_in) is False:
        print('file does not exist')
        return

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    sFilename_out= sFilename_output_in
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

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_discrete_in is not None:
        iFlag_discrete = iFlag_discrete_in
    else:
        iFlag_discrete = 0

    if iFlag_filter_in is not None:
        iFlag_filter = iFlag_filter_in
    else:
        iFlag_filter = 0

    if iFlag_thickness_in is not None:
        iFlag_thickness = iFlag_thickness_in
    else:
        iFlag_thickness = 0

    if sField_thickness_in is not None:
        sField_thickness = sField_thickness_in
    else:
        sField_thickness = ''

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if sField_color_in is not None:
        sField_color = sField_color_in
    else:
        sField_color = 'id'

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0


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

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'Spectral'

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title =1
    else:
        iFlag_title=0
        sTitle =  ''

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

    if iFlag_color == 1:
        pLayerdefn = pLayer.GetLayerDefn()
        nField = pLayerdefn.GetFieldCount()
        if nField == 0:
            iFlag_color = 0
            iFlag_field = 0
            iFlag_discrete = 0
        else:
            iFlag_field = 1
            sField0 = pLayerdefn.GetFieldDefn(0).name

        if sVariable_in is not None:
            sVariable = sVariable_in
        else:
            sVariable = sField0
    else:
        iFlag_field = 0
        iFlag_discrete = 0

    plt.rcParams["font.family"] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'
    fig = plt.figure( dpi=iDPI)
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    lID = 0
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aValue_field = list()

    if iFlag_color == 1:
        pLayerdefn = pLayer.GetLayerDefn()
        nField = pLayerdefn.GetFieldCount()
        if nField == 0:
            iFlag_color = 0
            iFlag_field = 0
            iFlag_discrete = 0
        else:
            iFlag_field = 1
            sField0 = pLayerdefn.GetFieldDefn(0).name

        if sVariable_in is not None:
            sVariable = sVariable_in
        else:
            sVariable = sField0
    else:
        iFlag_field = 0
        iFlag_discrete = 0

    nFeature = pLayer.GetFeatureCount()

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_field == 1:
            dValue = float(pFeature.GetField(sField_color))
            aValue_field.append(dValue)
        if sGeometry_type =='LINESTRING':
            aCoords_gcs =   get_geometry_coordinates(pGeometry_in)
            #aCoords_gcs = aCoords_gcs[:,0:2]
            dLon_max = float(np.max( [dLon_max, np.max(aCoords_gcs[:,0])] ))
            dLon_min = float(np.min( [dLon_min, np.min(aCoords_gcs[:,0])] ))
            dLat_max = float(np.max( [dLat_max, np.max(aCoords_gcs[:,1])] ))
            dLat_min = float(np.min( [dLat_min, np.min(aCoords_gcs[:,1])] ))

    if iFlag_field == 1:
        aValue_field = np.array(aValue_field)
        if iFlag_discrete ==1:
            #convert to integer
            aValue_field = np.array(aValue_field, dtype=int)
            #reorder it
            aValue_field = np.unique(aValue_field)
            aValue_field = np.sort(aValue_field)
            nValue_field = len(aValue_field) #get unique values from the aValue_field

        if iFlag_data_min == 1:  # min is provided
            dValue_min = dData_min  # np.min(aValue_field)
        else:
            aValue_field = aValue_field[aValue_field != dMissing_value]
            dValue_min = np.min(aValue_field)

        if iFlag_data_max == 1:  # max is provided
            dValue_max = dData_max  # np.max(aValue_field)
        else:
            aValue_field = aValue_field[aValue_field != dMissing_value]
            dValue_max = np.max(aValue_field)

        aValue_field = np.clip(aValue_field, a_min=dValue_min, a_max=dValue_max)

        print(dValue_min, dValue_max)

        # print(sVariable,dValue_min, dValue_max )
        if dValue_max == dValue_min:
            iFlag_same_value = 1
            return
        else:
            iFlag_same_value = 0
            pass

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = cpl.crs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                                central_latitude = 0.50*(dLat_max+dLat_min), globe=None)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84


    ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=pProjection_map ) #request.crs
    ax.set_global()


    if iFlag_discrete ==1:
        aIndex = np.linspace(0,1,nValue_field)
        prng = np.random.RandomState(1234567890)
        prng.shuffle(aIndex)
        #print(aIndex)
        colors = plt.colormaps[sColormap](aIndex)
        pCmap = ListedColormap(colors)
    else:
        pCmap = plt.colormaps[sColormap]

    if iFlag_thickness ==1:
        aValue_thickness =list()
        for pFeature in pLayer:
            dValue_thickness = pFeature.GetField(sField_thickness)
            aValue_thickness.append(dValue_thickness)

        aValue_thickness = np.array(aValue_thickness)
        dValue_thickness_max = np.max(aValue_thickness)
        dValue_thickness_min = np.min(aValue_thickness)

    iThickness_max = 2.5
    iThickness_min = 0.3
    #switch to collection in next development
    aPolyline = list()
    aColor = list()
    aThickness = list()
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in

    print(aExtent)
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    minx, miny, maxx, maxy = aExtent
    if iFlag_filter == 1:
        pLayer.SetSpatialFilterRect(minx, miny, maxx, maxy)

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_thickness ==1:
            dValue_thickness = pFeature.GetField(sField_thickness)
        if iFlag_color ==1:
            dValue_color = pFeature.GetField(sField_color)
        if sGeometry_type =='LINESTRING':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            aCoords_gcs = aCoords_gcs[:,0:2]
            nvertex = len(aCoords_gcs)
            if nvertex == 2 :
                dLon_label = 0.5 * (aCoords_gcs[0][0] + aCoords_gcs[1][0] )
                dLat_label = 0.5 * (aCoords_gcs[0][1] + aCoords_gcs[1][1] )
            else:
                lIndex_mid = int(nvertex/2)
                dLon_label = aCoords_gcs[lIndex_mid][0]
                dLat_label = aCoords_gcs[lIndex_mid][1]

            codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int )
            codes[0] = mpl.path.Path.MOVETO
            path = mpl.path.Path(aCoords_gcs, codes)
            x, y = zip(*path.vertices)

            if iFlag_thickness ==1:
                iThickness = remap( dValue_thickness, dValue_thickness_min, dValue_thickness_max, iThickness_min, iThickness_max )
            else:
                iThickness = 1.0

            aThickness.append(iThickness)

            if iFlag_color ==1:
                if iFlag_thickness ==1:
                    color_index = (dValue_thickness-dValue_thickness_min ) /(dValue_thickness_max - dValue_thickness_min )
                    color = pCmap(color_index)
                else:
                    iValue = int(dValue_color)
                    #find its index in the aValue_field array
                    iColor_index = np.where(aValue_field == iValue)[0][0]
                    color = pCmap(iColor_index)
            else:
                color = 'blue'

            aColor.append(color)
            aPolyline.append(list(zip(x, y)))
            lID = lID + 1

    pLC = LineCollection(aPolyline,  alpha=0.8, edgecolor=aColor,
                         facecolor='none', linewidths=aThickness, transform=pProjection_data)
    ax.add_collection(pLC)

    if iFlag_zebra ==1:
        ax.set_xticks(np.arange(minx, maxx+(maxx-minx)/11, (maxx-minx)/10))
        ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
        ax.set_axis_off()

    ax.coastlines(color='black', linewidth=1,resolution='10m')
    if iFlag_colorbar == 1:
        fig.canvas.draw()
        # Section 2
        ax_pos = ax.get_position() # get the original position
        #use this ax to set the colorbar ax position
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
                formatter = mpl.ticker.FuncFormatter(lambda x, pos: "{:.0f}".format(x))
                cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=pCmap,
                                           norm=mpl.colors.Normalize(
                                               dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
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

    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx+(maxx-minx)/9, (maxx-minx)/8),
                      ylocs=np.arange(miny, maxy+(maxy-miny)/9, (maxy-miny)/8))

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight': 'normal'}


    if iFlag_zebra ==1:
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size-2 )

            pass

    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass
        ax.set_title(sTitle)

    pDataset = pLayer = pFeature  = None
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    if sFilename_output_in is None:
        plt.show()
    else:
        #remove it if exists
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

