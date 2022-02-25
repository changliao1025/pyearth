
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
import cartopy.mpl.ticker as ticker
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

pProjection = ccrs.PlateCarree()
def map_raster_data(aImage_in, \
    aImage_extent, \
    sFilename_output_in,\
       iFlag_scientific_notation_colorbar_in=None,\
    sColormap_in = None,\
    iDPI_in = None,\
    dMissing_value_in=None,\
    dData_max_in = None, \
    dData_min_in = None,\
        sUnit_in=None):

    aImage_in = np.array(aImage_in)

    pShape = aImage_in.shape
    nrow, ncolumn = aImage_in.shape
    iSize_x = ncolumn
    iSize_y = nrow 
    sFilename_out= sFilename_output_in
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300
    
    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value= np.nanmin(aImage_in)

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
        sColormap =  'rainbow'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit =  ''    

    cmap = cm.get_cmap(sColormap)
   

    dummy_index = np.where(aImage_in > dData_max)
    aImage_in[dummy_index] = dData_max

    dummy_index = np.where(aImage_in < dData_min)
    aImage_in[dummy_index] = dData_min


    fig = plt.figure( dpi = iDPI  )
    #fig.set_figwidth( iSize_x )
    #fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.1, 0.65, 0.7], projection=pProjection )

    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)   

    rasterplot = ax.imshow(aImage_in, origin='upper', \
        extent=aImage_extent, \
        cmap = cmap, \
        transform=pProjection)   

    ax.coastlines(color='black', linewidth=1)
    
    ax.set_extent(aImage_extent)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k', 'rotation':90,'weight': 'normal'}
    ax_cb= fig.add_axes([0.75, 0.1, 0.02, 0.7])

    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \times 10^{{{}}}$'.format(a, b)
    
    if iFlag_scientific_notation_colorbar!=1:
        cb = plt.colorbar(rasterplot, cax = ax_cb, extend = 'max')
    else:
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        cb = plt.colorbar(rasterplot, cax = ax_cb, extend = 'max', format=formatter)
        
    cb.ax.get_yaxis().set_ticks_position('left')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel(sUnit, rotation=270)
    cb.ax.tick_params(labelsize=6) 

    plt.savefig(sFilename_out , bbox_inches='tight')
    #.show()

    plt.close('all')
    plt.clf()

    
    

