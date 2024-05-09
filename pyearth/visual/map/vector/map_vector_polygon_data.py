import os
import numpy as np
from osgeo import osr, gdal, ogr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter

#osr.UseExceptions()
#use agg and backend
#mpl.use('agg')
def map_vector_polygon_data(iFiletype_in,
                            sFilename_in,
                            sVariable_in=None,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFlag_color_in = None,
                            iFlag_colorbar_in=None,
                            iFlag_zebra_in=None,
                            iFlag_fill_in=None,
                            iFont_size_in=None,
                            iFlag_discrete_in=None,
                            sColormap_in=None,
                            sTitle_in=None,
                            iDPI_in=None,
                            iSize_x_in=None,
                            iSize_y_in=None,
                            dMissing_value_in=None,
                            dData_max_in=None,
                            dData_min_in=None,
                            sExtend_in=None,
                            sFont_in=None,
                            sUnit_in=None,
                            aLegend_in=None,
                            aExtent_in=None,
                            pProjection_map_in=None,
                            pProjection_data_in = None,
                            iFlag_debug = 0):
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

    if iFlag_fill_in is not None:
        iFlag_fill = iFlag_fill_in
    else:
        iFlag_fill = True

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

    if iFlag_discrete_in is not None:
        iFlag_discrete = iFlag_discrete_in
    else:
        iFlag_discrete = 0

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

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

    plt.rcParams["font.family"] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aValue = list()

    #get the field in the feature
    #pFeature = pLayer.GetNextFeature()

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

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_field == 1:
            dValue = float(pFeature.GetField(sVariable))
            aValue.append(dValue)
        if sGeometry_type == 'MULTIPOLYGON':
            for i in range(pGeometry_in.GetGeometryCount()):
                pPolygon = pGeometry_in.GetGeometryRef(i)
                aCoords_gcs = get_geometry_coordinates(pPolygon)
                dLon_max = np.max([dLon_max, np.max(aCoords_gcs[:, 0])])
                dLon_min = np.min([dLon_min, np.min(aCoords_gcs[:, 0])])
                dLat_max = np.max([dLat_max, np.max(aCoords_gcs[:, 1])])
                dLat_min = np.min([dLat_min, np.min(aCoords_gcs[:, 1])])
        else:
            if sGeometry_type == 'POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                dLon_max = np.max([dLon_max, np.max(aCoords_gcs[:, 0])])
                dLon_min = np.min([dLon_min, np.min(aCoords_gcs[:, 0])])
                dLat_max = np.max([dLat_max, np.max(aCoords_gcs[:, 1])])
                dLat_min = np.min([dLat_min, np.min(aCoords_gcs[:, 1])])

    if iFlag_field == 1:
        aValue = np.array(aValue)
        if iFlag_discrete ==1:
            #convert to integer
            aValue = np.array(aValue, dtype=int)
            #reorder it
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
        pProjection_map = ccrs.Orthographic(central_longitude=0.50*(
            dLon_max+dLon_min),  central_latitude=0.50*(dLat_max+dLat_min), globe=None)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        pProjection_data = pSRS_wgs84

        #pProjection_map = pSRS_wgs84
    #pProjection_map._threshold /= 1.0E6

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map)
    #ax.set_global()
    #use an advanced method for plotting
    if iFlag_discrete == 1:
        aIndex = np.linspace(0,1,nValue)
        prng = np.random.RandomState(1234567890)
        prng.shuffle(aIndex)
        #print(aIndex)
        colors = mpl.cm.get_cmap(sColormap)(aIndex)
        pCmap = ListedColormap(colors)
    else:
        pCmap = mpl.cm.get_cmap(sColormap)

    if aExtent_in is None:
        marginx = (dLon_max - dLon_min) / 20
        marginy = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx, dLon_max + marginx,
                   dLat_min - marginy, dLat_max + marginy]
    else:
        aExtent = aExtent_in

    print(aExtent)

    minx,  maxx, miny, maxy = aExtent
    #pLayer.SetSpatialFilterRect(minx, maxx, miny, maxy)
    ax.set_extent(aExtent, crs = pSRS_wgs84)
    ax.coastlines(linewidth=0.5, color='k', resolution='10m')
    if iFlag_debug ==1:
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            # get attribute
            if iFlag_field ==1:
                dValue = float(pFeature.GetField(sVariable))
                if dValue != dMissing_value:
                    if dValue > dValue_max:
                        dValue = dValue_max
                        #continue
                    if dValue < dValue_min:
                        dValue = dValue_min
                        #continue
                    if iFlag_discrete ==1:
                        #use unique value method to assign color
                        iValue = int(dValue)
                        #find its index in the aValue array
                        iColor_index = np.where(aValue == iValue)[0][0]
                        color = pCmap(iColor_index)
                    else:
                        iColor_index = int((dValue - dValue_min) /
                                       (dValue_max - dValue_min) * 255)
                        color = pCmap(iColor_index)
            else:
                color = 'blue'

            if sGeometry_type == 'MULTIPOLYGON':
                for i in range(pGeometry_in.GetGeometryCount()):
                    pPolygon = pGeometry_in.GetGeometryRef(i)
                    aCoords_gcs = get_geometry_coordinates(pPolygon)
                    aCoords_gcs = aCoords_gcs[:,0:2]
                    polygon = mpatches.Polygon(aCoords_gcs,  closed=True, fill=iFlag_fill,
                        alpha=0.8, edgecolor=None, facecolor=color,
                        transform= pProjection_data )
                    ax.add_patch(polygon)
                    #ax.plot(aCoords_gcs[:, 0], aCoords_gcs[:, 1], color=color,
                    #        transform=pSRS_wgs84)
            else:
                if sGeometry_type == 'POLYGON':
                    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    aCoords_gcs = aCoords_gcs[:,0:2]
                    polygon = mpatches.Polygon(aCoords_gcs, closed=True, fill=iFlag_fill,
                        alpha=0.8, edgecolor=None, facecolor=color,
                        transform= pProjection_data )
                    ax.add_patch(polygon)
                    #draw the polygon directly
                    #ax.plot(aCoords_gcs[:, 0], aCoords_gcs[:, 1], color=color,
                    #        transform=pSRS_wgs84)

    else:
        aPolygon = list()
        aColor_index=list()
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            # get attribute
            if iFlag_field == 1:
                dValue = float(pFeature.GetField(sVariable))
                if dValue != dMissing_value:
                    if dValue > dValue_max:
                        dValue = dValue_max
                        #continue
                    if dValue < dValue_min:
                        dValue = dValue_min
                        #continue
                    if iFlag_discrete ==1:
                        #use unique value method to assign color
                        iValue = int(dValue)
                        #find its index in the aValue array
                        iColor_index = np.where(aValue == iValue)[0][0]
                    else:
                        iColor_index = int((dValue - dValue_min) /
                                       (dValue_max - dValue_min) * 255)
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
            else:
                if sGeometry_type == 'POLYGON':
                    aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                    aCoords_gcs=aCoords_gcs[:,0:2]
                    aColor_index.append(iColor_index)
                    aPolygon.append(aCoords_gcs)


        aColor_index = np.array(aColor_index)
        #flatten the array as 1D
        aColor_index= aColor_index.flatten()
        aColor= pCmap(aColor_index)
        if iFlag_field == 1:
            aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
            if iFlag_fill == True:
                pPC = PatchCollection(aPatch, alpha=0.8,
                                      edgecolor=None,
                                      facecolor=aColor,
                                      linewidths=0.25,
                                      transform=pProjection_data)
            else:
                pPC = PatchCollection(aPatch, alpha=0.8,
                                      edgecolor=aColor,
                                      facecolor='none',
                                      linewidths=0.25,
                                      transform=pProjection_data)

        else:
            sColor = 'blue'
            aPatch = [Polygon(poly, closed=True, fill=iFlag_fill) for poly in aPolygon]
            if iFlag_fill == True:
                pPC = PatchCollection(aPatch, alpha=0.8,
                                      edgecolor=None,
                                      facecolor=sColor,
                                      linewidths=0.25,
                                      transform=pProjection_data)
            else:
                pPC = PatchCollection(aPatch, alpha=0.8,
                                      edgecolor=sColor,
                                      facecolor='none',
                                      linewidths=0.25,
                                      transform=pProjection_data)
        ax.add_collection(pPC)

    ax.set_extent(aExtent, crs = pSRS_wgs84)

    ax.set_xticks(np.arange(minx, maxx+(maxx-minx)/11, (maxx-minx)/10))
    ax.set_yticks(np.arange(miny, maxy+(maxy-miny)/11, (maxy-miny)/10))
    ax.set_axis_off()
    ax.set_title(sTitle)
    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=iFont_size-2)

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

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--',
                      xlocs=np.arange(minx, maxx, (maxx-minx)/4), ylocs=np.arange(miny, maxy, (maxy-miny)/4))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    if iFlag_zebra ==1:
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    pDataset = pLayer = pFeature = None
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
