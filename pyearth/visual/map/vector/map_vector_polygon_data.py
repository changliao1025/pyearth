
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.mpl.ticker as ticker
import matplotlib as mpl
from osgeo import  osr, gdal, ogr
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

pProjection = ccrs.PlateCarree()


class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format

def map_vector_polygon_data(iFiletype_in,\
                            sFilename_in, \
                            aImage_extent, \
                            sFilename_output_in,\
                            iFlag_scientific_notation_colorbar_in=None,\
                            iFlag_contour_in = None,\
                            sColormap_in = None,\
                            sTitle_in = None, \
                            iDPI_in = None,\
                            dMissing_value_in=None,\
                            dData_max_in = None, \
                            dData_min_in = None,\
                            sExtend_in =None,\
                            sUnit_in=None,\
                            aLegend_in = None):

    if iFiletype_in == 1: #geojson
        pDriver = ogr.GetDriverByName('GeoJSON')
    else:
        if iFiletype_in == 2: #shapefile
            pDriver = ogr.GetDriverByName('Esri Shapefile')

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)


    sFilename_out= sFilename_output_in
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

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'rainbow'

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

    cmap = cm.get_cmap(sColormap)



    fig = plt.figure( dpi = iDPI  )
    #fig.set_figwidth( iSize_x )
    #fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.1, 0.63, 0.7], projection=pProjection )

    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)
  

    ax.coastlines(color='black', linewidth=1)
    ax.set_title(sTitle)


    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = 0.06 + i * 0.04
            ax.text(0.03, dLocation, sText, \
                    verticalalignment='top', horizontalalignment='left',\
                    transform=ax.transAxes, \
                    color='black', fontsize=6)

            pass

    ax.set_extent(aImage_extent)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k', 'rotation':90,'weight': 'normal'}
    ax_cb= fig.add_axes([0.75, 0.1, 0.02, 0.7])





    plt.savefig(sFilename_out , bbox_inches='tight')
    #.show()

    plt.close('all')
    plt.clf()
