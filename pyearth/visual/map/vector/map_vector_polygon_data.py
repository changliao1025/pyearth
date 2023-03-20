import os
import numpy as np

import cartopy.crs as ccrs
import cartopy.mpl.ticker as ticker
import matplotlib as mpl
from shapely.wkt import loads
from osgeo import  osr, gdal, ogr

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

pProjection = ccrs.PlateCarree() #for latlon data only


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

def map_vector_polygon_data(iFiletype_in,
                            sFilename_in, 
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            sColormap_in = None,
                            sTitle_in = None, 
                            iDPI_in = None,
                            dMissing_value_in=None,
                            dData_max_in = None, 
                            dData_min_in = None,
                            sExtend_in =None,
                            sUnit_in=None,
                            aLegend_in = None,
                            aExtent_in = None,
                            pProjection_map_in=None):
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
    
    if iFiletype_in == 1: #geojson
        pDriver = ogr.GetDriverByName('GeoJSON')
    else:
        if iFiletype_in == 2: #shapefile
            pDriver = ogr.GetDriverByName('Esri Shapefile')

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)


    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

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
    iSize_x= 8 
    iSize_y= 8
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

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
        lID =0 
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                if dLon < dLon_min:
                    dLon_min = dLon
                if dLat > dLat_max:
                    dLat_max = dLat
                if dLat < dLat_min:
                    dLat_min = dLat

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),  central_latitude = 0.50*(dLat_max+dLat_min), globe=None)
   
    ax = fig.add_axes([0.1, 0.1, 0.63, 0.7], projection=pProjection_map )
    ax.set_global()
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        lID =0 
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25, \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)                   
    
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in
    
    ax.set_extent( aExtent )  
    

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

    

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k', 'rotation':90,'weight': 'normal'}
    sDirname = os.path.dirname(sFilename_output_in)

    pDataset = pLayer = pFeature  = None   
    if sFilename_output_in is None:
        plt.show()
    else:
        sFilename_out = os.path.join(sDirname, sFilename_output_in)
        plt.savefig(sFilename_out, bbox_inches='tight')   
        plt.close('all')
        plt.clf()
