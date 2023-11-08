import os
import numpy as np
from osgeo import osr, ogr
import matplotlib as mpl
import cartopy as cpl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
pProjection = cpl.crs.PlateCarree()  # for latlon data only


class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(
            self, useOffset=offset, useMathText=mathText)

    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom

    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format


def map_netcdf_file(sFilename_netcdf_in,
                    sVariable_in,
                    iFlag_unstructured_in=None,
                    sFolder_output_in=None,
                    iFlag_scientific_notation_colorbar_in=None,
                    iFont_size_in=None,
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
                    pBoundary_in=None):
    """
    Extract data from a NetCDF file based on a bounding box.

    Args:
        sFilename_netcdf_in (str): Path to the input NetCDF file.
        sFilename_netcdf_out (str): Path to the output NetCDF file.
        bbox (tuple): Bounding box as a tuple in the format (aLongitude_box_min, aLatitude_box_min, aLongitude_box_max, aLatitude_box_max).
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

    mpl.pyplot.rcParams["font.family"] = sFont

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
        ds = nc.Dataset(sFilename_netcdf_in, 'r')

        # get all the dimensions
        aDimension = ds.dimensions.keys()
        # get the time dimension by name
        aTime = ds.variables['time']
        nTime = len(aTime)

        aLatitude = ds.variables['lat'][:]
        aLongitude = ds.variables['lon'][:] - 180

        # read the variable
        aData = ds.variables[sVariable][:]
        # deal with the time

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
            dResolution_x = (dLongitude_max - dLongitude_min) / (ncolumn - 1)
            dResolution_y = (dLatitude_max - dLatitude_min) / (nrow - 1)

            # get max and min of the data

            if pBoundary_in is not None:

                aData_dump = list()
                for i in range(0, nrow):
                    for j in range(0, ncolumn):
                        dLongitude_cell_center = dLongitude_min + j * dResolution_x
                        dLatitude_cell_center = dLatitude_min + i * dResolution_y
                        # create
                        pPoint = ogr.Geometry(ogr.wkbPoint)
                        pPoint.AddPoint(dLongitude_cell_center,
                                        dLatitude_cell_center)
                        if pPoint.Within(pBoundary):
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

            # the bounding box should be applied to cell center
            for iStep in range(0, nTime):
                sStep = '{:03d}'.format(iStep)
                fig = mpl.pyplot.figure(dpi=iDPI)
                fig.set_figwidth(iSize_x)
                fig.set_figheight(iSize_y)
                ax = fig.add_axes([0.08, 0.1, 0.62, 0.7],
                                  projection=pProjection_map)
                ax.set_global()
                # redefine the filename and title
                sFilename_out = os.path.join(
                    sFolder_output, sVariable + '_' + sStep + '.png')

                for i in range(0, nrow):
                    for j in range(0, ncolumn):
                        dLongitude_cell_center = dLongitude_min + j * dResolution_x
                        dLatitude_cell_center = dLatitude_min + i * dResolution_y
                        # create
                        pPoint = ogr.Geometry(ogr.wkbPoint)
                        pPoint.AddPoint(dLongitude_cell_center,
                                        dLatitude_cell_center)
                        if pPoint.Within(pBoundary):
                            # get the value
                            dValue = aData[iStep, i, j]
                            # get its cell shape
                            aCoords_gcs = np.full((nVertex, 2), np.nan)

                            x1 = dLongitude_cell_center - dResolution_x/2
                            y1 = dLatitude_cell_center - dResolution_y/2
                            aCoords_gcs[0, :] = x1, y1

                            x1 = dLongitude_cell_center + dResolution_x/2
                            y1 = dLatitude_cell_center - dResolution_y/2
                            aCoords_gcs[1, :] = x1, y1

                            x1 = dLongitude_cell_center + dResolution_x/2
                            y1 = dLatitude_cell_center + dResolution_y/2
                            aCoords_gcs[2, :] = x1, y1

                            x1 = dLongitude_cell_center - dResolution_x/2
                            y1 = dLatitude_cell_center + dResolution_y/2
                            aCoords_gcs[3, :] = x1, y1

                            iColor_index = int(
                                (dValue - dValue_min) / (dValue_max - dValue_min) * 255)
                            # pick color from colormap
                            cmiColor_index = cmap(iColor_index)

                            polygon = mpl.patches.Polygon(aCoords_gcs[:, 0:2], closed=True, linewidth=0.25,
                                                       alpha=0.8, edgecolor=cmiColor_index, facecolor=cmiColor_index,
                                                       transform=cpl.crs.PlateCarree())
                            ax.add_patch(polygon)

                            pass
                        else:
                            pass

                # add extent
                if aExtent_in is None:
                    if pBoundary_in is None:
                        marginx = (dLongitude_max - dLongitude_min) / 20
                        marginy = (dLatitude_max - dLatitude_min) / 20
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

                if sFolder_output_in is None:
                    plt.show()
                else:

                    plt.savefig(sFilename_out, bbox_inches='tight')

                    plt.close('all')
                    plt.clf()

        print("Extraction successful!")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:


if __name__ == "__main__":
    sFilename_netcdf_in_path = "/compyfs/inputdata/lnd/dlnd7/mingpan/ming_daily_2019.nc"

    sVariable_in = "QOVER"

    sFolder_output = '/qfs/people/liao313/workspace/python/pyearth/figures/'

    aExtent = [-150.015625, -146.234375, 67.921875, 70.328125]
    dLongitude_left, dLongitude_right, dLatitude_bot, dLatitude_top = aExtent
    pRing = ogr.Geometry(ogr.wkbLinearRing)
    pRing.AddPoint(dLongitude_left, dLatitude_top)
    pRing.AddPoint(dLongitude_right, dLatitude_top)
    pRing.AddPoint(dLongitude_right, dLatitude_bot)
    pRing.AddPoint(dLongitude_left, dLatitude_bot)
    pRing.AddPoint(dLongitude_left, dLatitude_top)
    pBoundary = ogr.Geometry(ogr.wkbPolygon)
    pBoundary.AddGeometry(pRing)
    pBoundary_wkt = pBoundary.ExportToWkt()

    sTitle = 'Qver'
    sUnit = 'mm/day'

    map_netcdf_file(sFilename_netcdf_in_path, sVariable_in, sFolder_output_in=sFolder_output,
                    pBoundary_in=pBoundary_wkt, iFlag_scientific_notation_colorbar_in=1,
                    sTitle_in=sTitle, sUnit_in=sUnit)
