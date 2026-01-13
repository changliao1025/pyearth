import os
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt

plt.ion()
import cartopy as cpl
from osgeo import osr, gdal, ogr
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.system.define_global_variables import *
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.formatter import log_formatter
from pyearth.visual.formatter import OOMFormatter
from pyearth.visual.map.zebra_frame import zebra_frame
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates_batch

pProjection_default = cpl.crs.PlateCarree()


def map_raster_file(
    sFilename_in,
    sFilename_output_in=None,
    iFlag_scientific_notation_colorbar_in=None,
    iFlag_contour_in=None,
    iFlag_colorbar_in=None,
    iFlag_coastlines_in=None,
    sColormap_in=None,
    sTitle_in=None,
    iDPI_in=None,
    iFlag_zebra_in=None,
    dMissing_value_in=None,
    dData_max_in=None,
    dData_min_in=None,
    sExtend_in=None,
    sFormat_contour_in=None,
    sFont_in=None,
    sUnit_in=None,
    aExtent_in=None,
    pProjection_map_in=None,
    pProjection_data_in=None,
    aLabel_legend_in=None,
):

    pSRS_wgs84 = ccrs.PlateCarree()  # for latlon data only
    pSrs = osr.SpatialReference()
    # Initialize it with EPSG code for WGS84
    pSrs.ImportFromEPSG(4326)
    # Export to WKT
    pProjection_wgs84 = pSrs.ExportToWkt()
    pSRS_geodetic = ccrs.Geodetic()

    if os.path.exists(sFilename_in) == False:
        print("Error: file does not exist", sFilename_in)
        return

    # ask gdal to open the raster file
    pDataset = gdal.Open(sFilename_in)
    # get the number of bands
    iBand = pDataset.RasterCount
    # get the first band
    pBand = pDataset.GetRasterBand(1)
    # get the projection
    pProjection = pDataset.GetProjection()
    # get the geotransform
    aGeotransform = pDataset.GetGeoTransform()
    # get the extent
    dXmin = aGeotransform[0]
    dYmax = aGeotransform[3]
    dXmax = aGeotransform[0] + aGeotransform[1] * pDataset.RasterXSize
    dYmin = aGeotransform[3] + aGeotransform[5] * pDataset.RasterYSize
    aImage_extent = [dXmin, dXmax, dYmin, dYmax]
    # get the data of the first layer
    aImage_in = pBand.ReadAsArray()
    dNoData = pBand.GetNoDataValue()

    aImage_in = np.array(aImage_in).astype(float)

    pShape = aImage_in.shape
    nrow, ncolumn = aImage_in.shape
    iSize_x = ncolumn
    iSize_y = nrow
    sFilename_out = sFilename_output_in
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_contour_in is not None:
        iFlag_contour = iFlag_contour_in
    else:
        iFlag_contour = 0

    if iFlag_zebra_in is not None:
        iFlag_zebra = iFlag_zebra_in
    else:
        iFlag_zebra = 0

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        if dNoData is not None:
            dMissing_value = dNoData
        else:
            dMissing_value = np.nanmin(aImage_in)

    dummy_index = np.where(aImage_in == dMissing_value)
    if len(dummy_index[0]) > 0:
        aImage_in[dummy_index] = np.nan

    if dData_max_in is not None:
        dData_max = dData_max_in
    else:
        dData_max = np.nanmax(aImage_in)
        # print(dData_max)

    if dData_min_in is not None:
        dData_min = dData_min_in
    else:
        dData_min = np.nanmin(aImage_in)

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap = "rainbow"

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title = 1
    else:
        iFlag_title = 0
        sTitle = ""

    if sFormat_contour_in is not None:
        sFormat_contour = sFormat_contour_in
    else:
        sFormat_contour = "%1.1f"

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend = "max"

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit = ""

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams["font.family"] = "DeJavu Serif"
    plt.rcParams["font.serif"] = sFont
    plt.rcParams["mathtext.fontset"] = "dejavuserif"

    # cmap = mpl.cm.get_cmap(sColormap)
    cmap = plt.colormaps[sColormap]

    dummy_index = np.where(aImage_in > dData_max)
    aImage_in[dummy_index] = dData_max

    dummy_index = np.where(aImage_in < dData_min)
    aImage_in[dummy_index] = dData_min

    fig = plt.figure(dpi=iDPI)
    # fig.set_figwidth( iSize_x )
    # fig.set_figheight( iSize_y )

    pSpatial_refernce = osr.SpatialReference(pProjection)
    sCode = pSpatial_refernce.GetAttrValue("AUTHORITY", 1)

    if pProjection_data_in is not None:
        pProjection_data = pProjection_data_in
    else:
        if pSpatial_refernce.IsProjected():
            pProjection_data = ccrs.epsg(sCode)
        else:
            pProjection_data = pSRS_wgs84

    if aExtent_in is None:
        # check projection
        if pSpatial_refernce.IsProjected():
            # conver the extent to wgs84
            aX_new, aY_new = reproject_coordinates_batch(
                [aImage_extent[0], aImage_extent[1]],
                [aImage_extent[2], aImage_extent[3]],
                pSpatial_refernce.ExportToWkt(),
                pProjection_target_in=pProjection_wgs84,
            )
            aExtent_map = [aX_new[0], aX_new[1], aY_new[0], aY_new[1]]
        else:
            aExtent_map = aImage_extent
    else:
        aExtent_map = aExtent_in

    print(aExtent_map)
    minx, maxx, miny, maxy = aExtent_map  # these are in wgs84 projection

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = cpl.crs.Orthographic(
            central_longitude=0.50 * (maxx + minx),
            central_latitude=0.50 * (miny + maxy),
            globe=None,
        )

    ax = fig.add_axes([0.1, 0.1, 0.63, 0.7], projection=pProjection_map)

    # set a margin around the data
    # ax.set_xmargin(0.05)
    # ax.set_ymargin(0.10)
    ax.set_global()
    if iFlag_coastlines_in is not None:
        ax.coastlines(color="black", linewidth=1, resolution="10m")

    print(aImage_extent)
    ax.set_extent(aExtent_map, crs=pSRS_wgs84)
    rasterplot = ax.imshow(
        aImage_in,
        origin="upper",
        extent=aImage_extent,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
    )

    if iFlag_contour == 1:
        aPercentiles_in = np.arange(33, 67, 33)
        levels = cgpercentiles(aImage_in, aPercentiles_in, missing_value_in=-9999)
        contourplot = ax.contour(
            aImage_in,
            levels,
            colors="k",
            origin="upper",
            extent=aExtent_map,
            transform=pProjection_data,
            linewidths=0.5,
        )

        if iFlag_scientific_notation_colorbar == 1:
            ax.clabel(
                contourplot,
                contourplot.levels,
                inline=True,
                fmt=log_formatter,
                fontsize=7,
            )
        else:
            ax.clabel(
                contourplot,
                contourplot.levels,
                inline=True,
                fmt=sFormat_contour,
                fontsize=7,
            )

    ax.set_extent(aExtent_map, crs=pSRS_wgs84)
    # gridline
    gl = ax.gridlines(
        crs=cpl.crs.PlateCarree(),
        draw_labels=True,
        linewidth=1,
        color="gray",
        alpha=0.5,
        linestyle="--",
        xlocs=np.arange(minx, maxx + (maxx - minx) / 9, (maxx - minx) / 8),
        ylocs=np.arange(miny, maxy + (maxy - miny) / 9, (maxy - miny) / 8),
    )

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mpl.ticker.MaxNLocator(4)
    gl.ylocator = mpl.ticker.MaxNLocator(4)
    gl.xlabel_style = {"size": 10, "color": "k", "rotation": 0, "ha": "right"}
    gl.ylabel_style = {"size": 10, "color": "k", "rotation": 90, "weight": "normal"}

    rasterplot.set_clim(vmin=dData_min, vmax=dData_max)

    if iFlag_zebra == 1:
        ax.set_xticks(np.arange(minx, maxx + (maxx - minx) / 11, (maxx - minx) / 10))
        ax.set_yticks(np.arange(miny, maxy + (maxy - miny) / 11, (maxy - miny) / 10))
        ax.set_axis_off()

    if aLabel_legend_in is not None:
        # plot the first on the top
        sText = aLabel_legend_in[0]
        dLocation = 0.96
        ax.text(
            0.03,
            dLocation,
            sText,
            verticalalignment="top",
            horizontalalignment="left",
            transform=ax.transAxes,
            color="black",
            fontsize=10,
        )
        # plot the remaining on the bot
        nlegend = len(aLabel_legend_in)
        for i in range(1, nlegend, 1):
            sText = aLabel_legend_in[i]
            dLocation = nlegend * 0.06 - i * 0.05 - 0.03
            ax.text(
                0.03,
                dLocation,
                sText,
                verticalalignment="top",
                horizontalalignment="left",
                transform=ax.transAxes,
                color="black",
                fontsize=10,
            )

            pass

    if iFlag_title is None:
        ax.set_title(sTitle)
    else:
        if iFlag_title == 1:
            ax.set_title(sTitle)
        else:
            pass
        ax.set_title(sTitle)

    if iFlag_colorbar_in == 1:
        fig.canvas.draw()
        # Section 2
        ax_pos = ax.get_position()  # get the original position
        # use this ax to set the colorbar ax position
        ax_cb = fig.add_axes([ax_pos.x1 + 0.06, ax_pos.y0, 0.02, ax_pos.height])
        # ax_cb = fig.add_axes([0.75, 0.1, 0.02, 0.7])
        if iFlag_scientific_notation_colorbar == 1:
            formatter = OOMFormatter(fformat="%1.1e")
            cb = plt.colorbar(rasterplot, cax=ax_cb, extend=sExtend, format=formatter)
        else:
            formatter = OOMFormatter(fformat="%1.1f")
            cb = plt.colorbar(rasterplot, cax=ax_cb, extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position("right")
        cb.ax.get_yaxis().labelpad = 3
        cb.ax.set_ylabel(sUnit, rotation=90)
        cb.ax.get_yaxis().set_label_position("left")
        cb.ax.tick_params(labelsize=6)

    if iFlag_zebra == 1:
        ax.zebra_frame(crs=pSRS_wgs84, iFlag_outer_frame_in=1)

    ax.set_extent(aExtent_map, crs=pSRS_wgs84)

    if sFilename_output_in is None:
        plt.show()
        print("Finished plotting")
    else:
        # remove it if exists
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == ".png":
            plt.savefig(sFilename_out, bbox_inches="tight")
        else:
            if sExtension == ".pdf":
                plt.savefig(sFilename_out, bbox_inches="tight")
            else:
                plt.savefig(sFilename_out, bbox_inches="tight", format="ps")
        plt.close("all")
        plt.clf()

        print("Finish plotting raster map", sFilename_output_in)

    return
