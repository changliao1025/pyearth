import os
import sys
import numpy as np
from datetime import datetime
from osgeo import osr, ogr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import cartopy as cpl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.visual.formatter import OOMFormatter
pProjection = cpl.crs.PlateCarree()  # for latlon data only

def map_netcdf_file(sFilename_netcdf_in,
                    sVariable_in,
                    sFolder_output_in,
                    iFlag_monthly_in=None,
                    iFlag_daily_in=None,
                    iFlag_unstructured_in=None,                    
                    iFlag_scientific_notation_colorbar_in=None,
                    iFont_size_in=None,
                    sColormap_in=None,
                    sTitle_in=None,
                    iDPI_in=None,
                    iSize_x_in=None,
                    iSize_y_in=None,
                    dMissing_value_in=None,
                    dConvert_factor_in=None,
                    dData_max_in=None,
                    dData_min_in=None,
                    dResolution_x_in=None,
                    dResolution_y_in=None,
                    sExtend_in=None,
                    sFont_in=None,
                    sUnit_in=None,
                    aLegend_in=None,
                    aExtent_in=None,
                    pProjection_map_in=None,
                    pBoundary_in=None):
    """map a netcdf file using a variable

    Args:
        sFilename_netcdf_in (_type_): _description_
        sVariable_in (_type_): _description_
        sFolder_output_in (_type_): _description_
        iFlag_unstructured_in (_type_, optional): _description_. Defaults to None.
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        iFont_size_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        iSize_x_in (_type_, optional): _description_. Defaults to None.
        iSize_y_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dConvert_factor_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        dResolution_x_in (_type_, optional): _description_. Defaults to None.
        dResolution_y_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sFont_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
        pBoundary_in (_type_, optional): _description_. Defaults to None.

    Raises:
        ImportError: _description_
    """

    try:
        import netCDF4 as nc
    except ImportError as e:
        raise ImportError(
            "The package 'netCDF4' is required for this function to run.") from e

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 100

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

    if iFlag_monthly_in is not None:
        iFlag_monthly = iFlag_monthly_in
    else:
        iFlag_monthly = 0

    if iFlag_daily_in is not None:
        iFlag_daily = iFlag_daily_in
    else:
        iFlag_daily = 0

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

    if dData_min_in is not None:
        iFlag_data_min = 1
        dValue_min = dData_min_in
    else:
        iFlag_data_min = 0
        pass

    if dData_max_in is not None:
        iFlag_data_max = 1
        dValue_max = dData_max_in
    else:
        iFlag_data_max = 0
        pass

    if dConvert_factor_in is not None:
        dConvert_factor = dConvert_factor_in
    else:
        dConvert_factor = 1.0

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

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

    sVariable = sVariable_in

    if pBoundary_in is None:
        pBoundary = None
    else:
        # for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
        # https://gdal.org/api/python_gotchas.html
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)
        # pBoundary = pBoundary_in #only works for shapely geometry object
        pBoundary_box = pBoundary.GetEnvelope()

    cmap = mpl.cm.get_cmap(sColormap)

    # how can we define the boundary of the map?
    # if the boundary is not defined, we will use the boundary of the data
    # if the boundary is defined, we will use the boundary of the data

    # there are two ways to define the max and min of the data
    # the fist way is to use the max and min of the data,
    # this approach also have two sub-approaches: whole data, or data within the boundary
    # the second way is user defined max and min of the data

    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    try:
        pDatasets_in = nc.Dataset(sFilename_netcdf_in, 'r')
        pDatasets_in.set_auto_mask(False)

        #get the year information from the filename 
        sYear = sFilename_netcdf_in[-7:-3]
        iYear = int(sYear)
        # get all the dimensions
        
        for sKey, aValue in pDatasets_in.variables.items():    
            if sKey == 'time':
                aTime = aValue
                nTime = len(aTime)

            if sKey == 'lat':
                aLatitude = aValue[:]
                
            if sKey == 'lon':
                dummy = aValue[:]
                if np.min(dummy) < 0:
                    aLongitude = dummy
                    iFlag_cutoff = 0
                else:
                    aLongitude = dummy - 180
                    iFlag_cutoff = 1

            if sKey == sVariable:
                aData = aValue[:] * dConvert_factor
                
                pass

        
        
        # deal with the time
        print('Extracting data ...')

        if iFlag_unstructured_in is not None:
            # unstructured grid
            pass
        else:
            # regular grid
            nVertex = 4
            dLongitude_min = np.min(aLongitude)
            dLongitude_max = np.max(aLongitude)
            dLatitude_min = np.min(aLatitude)
            dLatitude_max = np.max(aLatitude)

            if pProjection_map_in is not None:
                pProjection_map = pProjection_map_in
            else:
                if pBoundary_in is None:
                    pProjection_map = cpl.crs.Orthographic(central_longitude=0.50*(dLongitude_min+dLongitude_max),
                                                        central_latitude=0.50*(dLatitude_min+dLatitude_max), globe=None)
                else:
                    # get the center of the boundary
                    centroid = pBoundary.Centroid()
                    # Extract coordinates
                    centroid_x = centroid.GetX()
                    centroid_y = centroid.GetY()
                    pProjection_map = cpl.crs.Orthographic(central_longitude=centroid_x,
                                                        central_latitude=centroid_y, globe=None)

                    pass

            ncolumn = len(aLongitude)
            nrow = len(aLatitude)
            if iFlag_cutoff == 1:
                #shift the data
                shift = int(ncolumn/2)
                aData = np.roll(aData, shift, axis=2)                   
                pass
            
            if dResolution_x_in is not None:
                dResolution_x = dResolution_x_in
            else:
                dResolution_x = (dLongitude_max - dLongitude_min) / (ncolumn - 1)
            
            if dResolution_y_in is not None:
                dResolution_y = dResolution_y_in
            else:
                dResolution_y = (dLatitude_max - dLatitude_min) / (nrow - 1)

            # get max and min of the data
            print(dResolution_x, dResolution_y)
            # Pre-calculate half resolutions
            half_res_x = dResolution_x / 2
            half_res_y = dResolution_y / 2


            #predefine a lat-lon meshgrid using the center of the cell
            aLongitude_center = np.arange(ncolumn) * dResolution_x + dLongitude_min + half_res_x
            aLatitude_center = np.arange(nrow) * dResolution_y + dLatitude_min + half_res_y

            aLongitude_grid, aLatitude_grid = np.meshgrid(aLongitude_center, aLatitude_center)
            # Create a mask for cells within the boundary
            aMask = np.ones((nrow, ncolumn), dtype=bool)
            if pBoundary_in is not None:
                for i in range(nrow):
                    for j in range(ncolumn):
                        pPoint = ogr.Geometry(ogr.wkbPoint)
                        pPoint.AddPoint(aLongitude_grid[i, j], aLatitude_grid[i, j])
                        if not pPoint.Within(pBoundary):
                            aMask[i, j] = False

            #find min and max of the data
            if iFlag_data_max == 1 and iFlag_data_min == 1:
                aData_dump = list()
                for i in range(0, nrow):
                    for j in range(0, ncolumn):
                        if not aMask[i, j]:
                            continue
                        aData_dump.append(aData[:, i, j])
                        pass
                aData_dump = np.array(aData_dump)
                print(np.max(aData_dump))
                pass
            else:
                if pBoundary_in is not None:
                    aData_dump = list()
                    for i in range(0, nrow):
                        for j in range(0, ncolumn):
                            if not aMask[i, j]:
                                continue

                            aData_dump.append(aData[:, i, j])
                            pass

                    aData_dump = np.array(aData_dump)
                    if iFlag_data_max == 0:
                        dValue_max = np.max(aData_dump)
                    if iFlag_data_min == 0:
                        dValue_min = np.min(aData_dump)
                else:
                    if iFlag_data_max == 0:
                        dValue_max = np.max(aData)

                    if iFlag_data_min == 0:
                        dValue_min = np.min(aData)

            print(dValue_min, dValue_max)
            #flush the print buffer 
            
            sys.stdout.flush()
            
            # the bounding box should be applied to cell center
            
            #use the first day as index 
            pDate0 = datetime(iYear, 1, 1)
            if iFlag_monthly == 1:
                for iMonth in range(1,13):
                    sMonth = '{:02d}'.format(iMonth)
                    fig = plt.figure(dpi=iDPI)
                    fig.set_figwidth(iSize_x)
                    fig.set_figheight(iSize_y)
                    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7],
                                      projection=pProjection_map)
                    ax.set_global()
                    # redefine the filename and title
                    sFilename_out = os.path.join(
                        sFolder_output_in, sVariable + '_' + sMonth + '.png')

                    #get the start and end index 
                    pDate_start = datetime(iYear, iMonth, 1)
                    if iMonth == 12:
                        pDate_end = datetime(iYear+1, 1, 1)
                    else:
                        pDate_end = datetime(iYear, iMonth+1, 1)

                    iStart = (pDate_start - pDate0).days
                    iEnd = (pDate_end - pDate0).days
                    
                    aPolygon = []
                    aColor = []
                    aData_dummy = []
                    for i in range(0, nrow):
                        for j in range(0, ncolumn):
                            if not aMask[i, j]:
                                continue

                            dLongitude_cell_center = aLongitude_grid[i, j]
                            dLatitude_cell_center = aLatitude_grid[i, j]
                            data_dummy = aData[iStart:iEnd, i, j]
                            dValue = np.mean(data_dummy)
                            aData_dummy.append(dValue)

                            # get its cell shape
                            aCoords_gcs = np.full((nVertex, 2), np.nan)

                            # Define coordinates for each vertex
                            coords = [(dLongitude_cell_center - half_res_x, dLatitude_cell_center - half_res_y),
                                      (dLongitude_cell_center + half_res_x, dLatitude_cell_center - half_res_y),
                                      (dLongitude_cell_center + half_res_x, dLatitude_cell_center + half_res_y),
                                      (dLongitude_cell_center - half_res_x, dLatitude_cell_center + half_res_y)]

                            # Fill aCoords_gcs using a loop
                            for idx, (x, y) in enumerate(coords):
                                aCoords_gcs[idx, :] = x, y

                            #dealing with missing value
                            if dValue > dValue_min: 
                                iColor_index = np.clip(((dValue - dValue_min) / (dValue_max - dValue_min) * 255).astype(int), 0, 255)
                                # pick color from colormap       
                                aColor.append(cmap(iColor_index))
                                aPolygon.append(aCoords_gcs[:, 0:2])

                    print(np.max(aData_dummy))
                    # add extent
                    # Create a PatchCollection with the polygons
                    aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
                    pPC = PatchCollection(aPatch, cmap=cmap, alpha=0.8, edgecolor=None, 
                                          facecolor=aColor, linewidths=0.25, 
                                          transform=cpl.crs.PlateCarree())
                    ax.add_collection(pPC)
                    if aExtent_in is None:
                        if pBoundary_in is None:
                            marginx = (dLongitude_max - dLongitude_min) / 20
                            marginy = (dLatitude_max - dLatitude_min) / 20
                            if dLongitude_min - marginx  < -180 or dLongitude_max + marginx > 180:
                                marginx = 0                       
                            if dLatitude_min - marginy < -90 or dLatitude_max + marginy > 90:
                                marginy = 0                       

                            aExtent = [dLongitude_min - marginx, dLongitude_max +
                                       marginx, dLatitude_min - marginy, dLatitude_max + marginy]
                        else:

                            marginx = (pBoundary_box[1] - pBoundary_box[0]) / 20
                            marginy = (pBoundary_box[3] - pBoundary_box[2]) / 20
                            aExtent = [pBoundary_box[0] - marginx, pBoundary_box[1] +
                                       marginx, pBoundary_box[2] - marginy, pBoundary_box[3] + marginy]
                    else:
                        aExtent = aExtent_in

                    ax.set_extent(aExtent)
                    ax.coastlines(color='black', linewidth=1)
                    ax.set_title(sTitle)
                    ax_cb = fig.add_axes([0.75, 0.15, 0.02, 0.6])

                    if iFlag_scientific_notation_colorbar == 1:
                        formatter = OOMFormatter(fformat="%1.1e")
                        cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                                       cmap=cmap,
                                                       norm=mpl.colors.Normalize(
                                                           dValue_min, dValue_max),  # vmax and vmin
                                                       extend=sExtend, format=formatter)
                    else:
                        formatter = OOMFormatter(fformat="%1.2f")
                        cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                                       cmap=cmap,
                                                       norm=mpl.colors.Normalize(
                                                           dValue_min, dValue_max),  # vmax and vmin
                                                       extend=sExtend, format=formatter)

                    cb.ax.get_yaxis().set_ticks_position('right')
                    cb.ax.get_yaxis().labelpad = 5
                    cb.ax.set_ylabel(sUnit, rotation=90, fontsize=iFont_size-2)
                    cb.ax.get_yaxis().set_label_position('left')
                    cb.ax.tick_params(labelsize=iFont_size-2)

                    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER

                    gl.xlabel_style = {'size': 10, 'color': 'k',
                                       'rotation': 0, 'ha': 'right'}
                    gl.ylabel_style = {'size': 10, 'color': 'k',
                                       'rotation': 90, 'weight': 'normal'}



                    plt.savefig(sFilename_out, bbox_inches='tight')
                    plt.close('all')
                    plt.clf()
                    print(sFilename_out)
                    pass
                pass
            else:
                pass
            if iFlag_daily == 1:
                # Pre-calculate the color indices
                color_indices = np.clip(((aData - dValue_min) / (dValue_max - dValue_min) * 255).astype(int), 0, 255)
            

                for iStep in range(0, nTime):
                    sStep = '{:03d}'.format(iStep+1)
                    fig = plt.figure(dpi=iDPI)
                    fig.set_figwidth(iSize_x)
                    fig.set_figheight(iSize_y)
                    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7],
                                      projection=pProjection_map)
                    ax.set_global()
                    # redefine the filename and title
                    sFilename_out = os.path.join(
                        sFolder_output_in, sVariable + '_' + sStep + '.png')

                    aPolygon = []
                    aColor = []
                    for i in range(0, nrow):
                        for j in range(0, ncolumn):
                            if not aMask[i, j]:
                                continue

                            dLongitude_cell_center = aLongitude_grid[i, j]
                            dLatitude_cell_center = aLatitude_grid[i, j]
                            dValue = aData[iStep, i, j]

                            # get its cell shape
                            aCoords_gcs = np.full((nVertex, 2), np.nan)

                            # Define coordinates for each vertex
                            coords = [(dLongitude_cell_center - half_res_x, dLatitude_cell_center - half_res_y),
                                      (dLongitude_cell_center + half_res_x, dLatitude_cell_center - half_res_y),
                                      (dLongitude_cell_center + half_res_x, dLatitude_cell_center + half_res_y),
                                      (dLongitude_cell_center - half_res_x, dLatitude_cell_center + half_res_y)]

                            # Fill aCoords_gcs using a loop
                            for idx, (x, y) in enumerate(coords):
                                aCoords_gcs[idx, :] = x, y

                            #dealing with missing value
                            if dValue > dValue_min:
                                iColor_index = color_indices[iStep, i, j]
                                # pick color from colormap       
                                aColor.append(cmap(iColor_index))
                                aPolygon.append(aCoords_gcs[:, 0:2])

                    # add extent
                    # Create a PatchCollection with the polygons
                    aPatch = [Polygon(poly, closed=True) for poly in aPolygon]
                    pPC = PatchCollection(aPatch, cmap=cmap, alpha=0.8, edgecolor=None, 
                                          facecolor=aColor, linewidths=0.25, 
                                          transform=cpl.crs.PlateCarree())
                    ax.add_collection(pPC)
                    if aExtent_in is None:
                        if pBoundary_in is None:
                            marginx = (dLongitude_max - dLongitude_min) / 20
                            marginy = (dLatitude_max - dLatitude_min) / 20
                            if dLongitude_min - marginx  < -180 or dLongitude_max + marginx > 180:
                                marginx = 0                       
                            if dLatitude_min - marginy < -90 or dLatitude_max + marginy > 90:
                                marginy = 0                       

                            aExtent = [dLongitude_min - marginx, dLongitude_max +
                                       marginx, dLatitude_min - marginy, dLatitude_max + marginy]
                        else:

                            marginx = (pBoundary_box[1] - pBoundary_box[0]) / 20
                            marginy = (pBoundary_box[3] - pBoundary_box[2]) / 20
                            aExtent = [pBoundary_box[0] - marginx, pBoundary_box[1] +
                                       marginx, pBoundary_box[2] - marginy, pBoundary_box[3] + marginy]
                    else:
                        aExtent = aExtent_in

                    ax.set_extent(aExtent)
                    ax.coastlines(color='black', linewidth=1)
                    ax.set_title(sTitle)
                    ax_cb = fig.add_axes([0.75, 0.15, 0.02, 0.6])

                    if iFlag_scientific_notation_colorbar == 1:
                        formatter = OOMFormatter(fformat="%1.1e")
                        cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                                       cmap=cmap,
                                                       norm=mpl.colors.Normalize(
                                                           dValue_min, dValue_max),  # vmax and vmin
                                                       extend=sExtend, format=formatter)
                    else:
                        formatter = OOMFormatter(fformat="%1.2f")
                        cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                                       cmap=cmap,
                                                       norm=mpl.colors.Normalize(
                                                           dValue_min, dValue_max),  # vmax and vmin
                                                       extend=sExtend, format=formatter)

                    cb.ax.get_yaxis().set_ticks_position('right')
                    cb.ax.get_yaxis().labelpad = 5
                    cb.ax.set_ylabel(sUnit, rotation=90, fontsize=iFont_size-2)
                    cb.ax.get_yaxis().set_label_position('left')
                    cb.ax.tick_params(labelsize=iFont_size-2)

                    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
                    gl.xformatter = LONGITUDE_FORMATTER
                    gl.yformatter = LATITUDE_FORMATTER

                    gl.xlabel_style = {'size': 10, 'color': 'k',
                                       'rotation': 0, 'ha': 'right'}
                    gl.ylabel_style = {'size': 10, 'color': 'k',
                                       'rotation': 90, 'weight': 'normal'}



                    plt.savefig(sFilename_out, bbox_inches='tight')
                    plt.close('all')
                    plt.clf()
                    print(sFilename_out)

        print("Extraction successful!")
    except Exception as e:
        print(f"An error occurred: {e}")


