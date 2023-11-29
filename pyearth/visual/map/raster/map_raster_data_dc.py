
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy as cpl

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.visual.formatter import log_formatter
pProjection = cpl.crs.PlateCarree()



def map_raster_data_dc(aImage_in,
                       aImage_extent,
                       sFilename_output_in,
                       iFlag_scientific_notation_colorbar_in=None,
                       iFlag_contour_in=None,
                       sTitle_in=None,
                       aInterval_in=None,
                       aColor_in=None,
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

    if sFormat_contour_in is not None:
        sFormat_contour = sFormat_contour_in
    else:
        sFormat_contour = '%1.1f'

    if aInterval_in is not None:
        iFlag_interval = 1
        aInterval = aInterval_in
        nInterval = len(aInterval)
    else:
        aPercentiles_in = np.arange(33, 67, 33)
        aInterval = cgpercentiles(
            aImage_in, aPercentiles_in, missing_value_in=dMissing_value_in)
        iFlag_interval = 0

    if aColor_in is not None:
        iFlag_colar = 1
        aColor = aColor_in
        ncolor = len(aColor)
        if iFlag_interval == 1:
            ncolor = nInterval + 1
        else:
            pass
    else:
        iFlag_colar = 0
        if iFlag_interval == 1:
            ncolor = nInterval + 1
            aColor = create_diverge_rgb_color_hex(ncolor)
        else:

            pass

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

    if dData_min_in is not None:
        dData_min = dData_min_in
    else:
        dData_min = np.nanmin(aImage_in)

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

    dummy_index = np.where(aImage_in > dData_max)
    aImage_in[dummy_index] = dData_max

    dummy_index = np.where(aImage_in < dData_min)
    aImage_in[dummy_index] = dData_min

    fig = plt.figure(dpi=iDPI)

    ax = fig.add_axes([0.1, 0.1, 0.63, 0.7], projection=pProjection)

    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    if iFlag_contour == 1:
        aPercentiles_in = np.arange(33, 67, 33)
        levels = cgpercentiles(aImage_in, aPercentiles_in,
                               missing_value_in=dMissing_value)
        contourplot = ax.contour(aImage_in, levels, colors='k', origin='upper',
                                 extent=aImage_extent, transform=pProjection, linewidths=0.5)

        if iFlag_scientific_notation_colorbar == 1:
            ax.clabel(contourplot, contourplot.levels,
                      inline=True, fmt=log_formatter, fontsize=7)
        else:
            ax.clabel(contourplot, contourplot.levels,
                      inline=True, fmt=sFormat_contour, fontsize=7)

    aPseudo_image = np.full((nrow, ncolumn), fill_value=np.nan)

    dummy_index = np.where(aImage_in <= aInterval[0])
    if len(dummy_index) != 0:
        aPseudo_image[dummy_index] = 0

    for i in range(0, nInterval-1, 1):
        dIntervel0 = aInterval[i]
        dIntervel1 = aInterval[i+1]
        dummy_index = np.where((aImage_in > dIntervel0)
                               & (aImage_in <= dIntervel1))
        if len(dummy_index) != 0:
            aPseudo_image[dummy_index] = i+1

    dummy_index = np.where(aImage_in > aInterval[nInterval-1])
    if len(dummy_index) != 0:
        aPseudo_image[dummy_index] = nInterval

    cmap = (mpl.colors.ListedColormap(aColor[1:ncolor-1])
            .with_extremes(over=aColor[ncolor-1], under=aColor[0]))

    uni = np.unique(aPseudo_image[np.where(np.isfinite(aPseudo_image))])
    uni_index = np.array(uni).astype(int)
    dummy_color = list()
    for i in uni_index:
        dummy_color.append(aColor[i])
    cmap0 = mpl.colors.ListedColormap(dummy_color)

    dcrasterplot = ax.imshow(aPseudo_image, origin='upper',
                             extent=aImage_extent,
                             cmap=cmap0,
                             transform=pProjection)

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

    ax.coastlines(color='black', linewidth=1)
    ax.set_title(sTitle)
    ax.set_extent(aImage_extent)

    gl = ax.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    ax_cb = fig.add_axes([0.75, 0.1, 0.02, 0.7])
    bounds = np.linspace(0, nInterval-1, nInterval, endpoint=True)

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), cax=ax_cb,
                      extend=sExtend, extendfrac='auto',
                      ticks=bounds)

    cb.ax.get_yaxis().set_ticks_position('right')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel(sUnit, rotation=270)
    cb.ax.tick_params(labelsize=6)
    aLabel = list()
    for i in range(nInterval):
        if iFlag_scientific_notation_colorbar == 1:
            sLabel = '{:.1e}'.format(aInterval[i])
        else:
            sLabel = '{:.1f}'.format(aInterval[i])
        aLabel.append(sLabel)

    bounds0 = np.linspace(0, nInterval-1, nInterval)
    bounds = bounds0[0:nInterval]
    cb.ax.set_yticks(bounds)
    cb.ax.set_yticklabels(aLabel)
    # cb.ax.set_yticklabels([aLabel[int(i)] for i in bounds]) # add the labels

    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close('all')
    plt.clf()
