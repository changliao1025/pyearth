
import pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

pProjection = ccrs.PlateCarree()
def map_raster_data(aImage_in, \
    aImage_extent, \
    sFilename_output_in,\
    sColormap_in = None,\
    iDPI_in = None,\
    dMissing_value_in=None,\
    dData_max_in = None, \
    dData_min_in = None):

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

    cmap = cm.get_cmap(sColormap)
   

    dummy_index = np.where(aImage_in > dData_max)
    aImage_in[dummy_index] = dData_max

    dummy_index = np.where(aImage_in < dData_min)
    aImage_in[dummy_index] = dData_min


    fig = plt.figure( dpi = iDPI  )
    #fig.set_figwidth( iSize_x )
    #fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.7], projection=pProjection )


    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)   

    ax.imshow(aImage_in, origin='upper', \
        extent=aImage_extent, \
        cmap = cmap, \
        transform=pProjection)   

    ax.coastlines(color='black', linewidth=1)
    #ax.set_global()
    #ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    #ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())    
    #lon_formatter = LongitudeFormatter(number_format='.1f',
    #                                   degree_symbol='',
    #                                   dateline_direction_label=True)
    #lat_formatter = LatitudeFormatter(number_format='.1f',
    #                                  degree_symbol='')
    #ax.xaxis.set_major_formatter(lon_formatter)
    #ax.yaxis.set_major_formatter(lat_formatter)
    ax.set_extent(aImage_extent)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

    plt.savefig(sFilename_out , bbox_inches='tight')
    #.show()

    plt.close('all')
    plt.clf()

    
    

