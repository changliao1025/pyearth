
import numpy as np


import matplotlib as mpl
import cartopy as cpl

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.toolbox.data.cgpercentiles import cgpercentiles


pProjection = cpl.crs.PlateCarree()


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


def map_vector_point_data(iFiletype_in,
                          sFilename_in,
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
                          sUnit_in=None,
                          aLegend_in=None):

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

    fig = mpl.pyplot.figure(dpi=iDPI)
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
                      inline=True, fmt=fmt0, fontsize=4)
        else:

            ax.clabel(contourplot, contourplot.levels,
                      inline=True, fmt=fmt1, fontsize=4)

    ax.coastlines(color='black', linewidth=1)
    ax.set_title(sTitle)

    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = 0.06 + i * 0.04
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=6)

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
        cb = mpl.pyplot.colorbar(rasterplot, cax=ax_cb,
                          extend=sExtend, format=formatter)
    else:
        formatter = OOMFormatter(fformat="%1.1f")
        cb = mpl.pyplot.colorbar(rasterplot, cax=ax_cb,
                          extend=sExtend, format=formatter)

    cb.ax.get_yaxis().set_ticks_position('right')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel(sUnit, rotation=270)
    cb.ax.tick_params(labelsize=6)

    mpl.pyplot.savefig(sFilename_out, bbox_inches='tight')
    # .show()

    mpl.pyplot.close('all')
    mpl.pyplot.clf()
