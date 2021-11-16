
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs

pProjection = ccrs.PlateCarree()
def map_raster_data(aImage_in, aImage_extent, \
    sFilename_output_in,\
    sColormap_in = None,\
    iDPI_in = None,\
    missing_value_in=None,\
    data_max_in = None, \
    data_min_in = None):

    aImage_in = np.array(aImage_in)

    pShape = aImage_in.shape
    nrow, ncolumn = aImage_in.shape
    iSize_x = ncolumn * 2
    iSize_y = nrow * 2

    sFilename_out= sFilename_output_in
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if data_max_in is not None:
        data_max = data_max_in
    else:
        data_max = np.nanmax(aImage_in)

    if data_min_in is not None:
        data_min = data_min_in
    else:
        data_min = np.nanmin(aImage_in)
    
    if missing_value_in is not None:
        missing_value = missing_value_in
    else:
        missing_value = np.nanmin(aImage_in)

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'rainbow'

    cmap = cm.get_cmap(sColormap)

    dummy_index = np.where(aImage_in > data_max)
    aImage_in[dummy_index] = data_max

    dummy_index = np.where(aImage_in < data_min)
    aImage_in[dummy_index] = data_min


    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4], projection=pProjection )


    # set a margin around the data
    ax.set_xmargin(0.05)
    ax.set_ymargin(0.10)

    ax.imshow(aImage_in, origin='upper', \
        extent=aImage_extent, \
        cmap = cmap, \
        transform=ccrs.PlateCarree())

    ax.coastlines(resolution='50m', color='black', linewidth=1)

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()

    
    

