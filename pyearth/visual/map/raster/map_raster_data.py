
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy as cpl

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.formatter import log_formatter
from pyearth.visual.formatter import OOMFormatter

pProjection = cpl.crs.PlateCarree()

def map_raster_data(aImage_in,
                    aImage_extent,
                    sFilename_output_in,
                    iFlag_scientific_notation_colorbar_in=None,
                    iFlag_contour_in=None,
                    sColormap_in=None,
                    sTitle_in=None,
                    iDPI_in=None,
                    dMissing_value_in=None,
                    dData_max_in=None,
                    dData_min_in=None,
                    sExtend_in=None,
                    sFormat_contour_in=None,
                    sUnit_in=None,
                    aLabel_legend_in=None):

    aImage_in = np.array(aImage_in)

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

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = np.nanmin(aImage_in)

    dummy_index = np.where(aImage_in == dMissing_value)
    aImage_in[dummy_index] = np.nan

    if dData_max_in is not None:
        dData_max = dData_max_in
    else:
        dData_max = np.nanmax(aImage_in)
        print(dData_max)

    if dData_min_in is not None:
        dData_min = dData_min_in
    else:
        dData_min = np.nanmin(aImage_in)

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

    if sFormat_contour_in is not None:
        sFormat_contour = sFormat_contour_in
    else:
        sFormat_contour = '%1.1f'

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend = 'max'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit = ''

    cmap = mpl.cm.get_cmap(sColormap)

    dummy_index = np.where(aImage_in > dData_max)
    aImage_in[dummy_index] = dData_max

    dummy_index = np.where(aImage_in < dData_min)
    aImage_in[dummy_index] = dData_min

    fig = plt.figure(dpi=iDPI)
    # fig.set_figwidth( iSize_x )
    # fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.1, 0.63, 0.7], projection=pProjection)

    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    rasterplot = ax.imshow(aImage_in, origin='upper',
                           extent=aImage_extent,
                           cmap=cmap,
                           transform=pProjection)

    if iFlag_contour == 1:
        aPercentiles_in = np.arange(33, 67, 33)
        levels = cgpercentiles(
            aImage_in, aPercentiles_in, missing_value_in=-9999)
        contourplot = ax.contour(aImage_in, levels, colors='k', origin='upper',
                                 extent=aImage_extent, transform=pProjection, linewidths=0.5)

        if iFlag_scientific_notation_colorbar == 1:
            ax.clabel(contourplot, contourplot.levels,
                      inline=True, fmt=log_formatter, fontsize=7)
        else:
            ax.clabel(contourplot, contourplot.levels,
                      inline=True, fmt=sFormat_contour, fontsize=7)

    ax.coastlines(color='black', linewidth=1)
    ax.set_title(sTitle)

    if aLabel_legend_in is not None:
        # plot the first on the top
        sText = aLabel_legend_in[0]
        dLocation = 0.96
        ax.text(0.03, dLocation, sText,
                verticalalignment='top', horizontalalignment='left',
                transform=ax.transAxes,
                color='black', fontsize=10)
        # plot the remaining on the bot
        nlegend = len(aLabel_legend_in)
        for i in range(1, nlegend, 1):
            sText = aLabel_legend_in[i]
            dLocation = nlegend * 0.06 - i * 0.05 - 0.03
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=10)

            pass

    ax.set_extent(aImage_extent)

    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}
    ax_cb = fig.add_axes([0.75, 0.1, 0.02, 0.7])

    rasterplot.set_clim(vmin=dData_min, vmax=dData_max)

    if iFlag_scientific_notation_colorbar == 1:
        formatter = OOMFormatter(fformat="%1.1e")
        cb = plt.colorbar(rasterplot, cax=ax_cb,
                          extend=sExtend, format=formatter)
    else:
        formatter = OOMFormatter(fformat="%1.1f")
        cb = plt.colorbar(rasterplot, cax=ax_cb,
                          extend=sExtend, format=formatter)

    cb.ax.get_yaxis().set_ticks_position('right')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel(sUnit, rotation=270)
    cb.ax.tick_params(labelsize=6)

    plt.savefig(sFilename_out, bbox_inches='tight')
    # .show()

    plt.close('all')
    plt.clf()
