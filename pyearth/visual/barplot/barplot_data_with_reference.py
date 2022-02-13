import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker


from pyearth.system.define_global_variables import *

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def barplot_data_with_reference(aData_in,\
                 aLabel_x_in,\
                 aLabel_y_in,\
                 sFilename_out,\
                     aReference_in,\
                 iDPI_in = None,\
                 ncolumn_in = None,\
                 iSize_x_in = None, \
                 iSize_y_in = None, \
                 dMax_y_in = None, \
                 dMin_y_in = None, \
                 dSpace_y_in = None,\
                 aMarker_in = None,\
                 aColor_in = None,\
                 aHatch_in = None,\
                     sLabel_info_in = None,\
                 sLabel_y_in = None, \
                 aLabel_legend_in = None,\
                     aLinestyle_in = None,\
                 aLocation_legend_in =None,\
                 sFormat_y_in =None,\
                 sLocation_legend_in=None,\
                 sTitle_in = None):

    aData_in = np.array(aData_in)
    pShape = aData_in.shape

    nData = pShape[0] #len(aLabel_x_in)

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 12

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 9

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = np.nanmax(aData_in) * 1.0

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aData_in) #if it has negative value, change here

    if (dMax_y <= dMin_y ):
        return

    if sLabel_y_in is not None:
        sLabel_y = sLabel_y_in
    else:
        sLabel_y = ''

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    if sFormat_y_in is not None:
        iFlag_format_y = 1
        sFormat_y = sFormat_y_in
    else:
        iFlag_format_y = 0

    if sLocation_legend_in is not None:
        sLocation_legend = sLocation_legend_in
    else:
        sLocation_legend = "upper right"

    if aLocation_legend_in is not None:
        aLocation_legend = aLocation_legend_in
    else:
        aLocation_legend=(1.0,1.0)

    if ncolumn_in is not None:
        ncolumn = ncolumn_in
    else:
        ncolumn = 1

    if aColor_in is not None:
        aColor = aColor_in
    else:
        if(nData>=3):
            aColor= create_diverge_rgb_color_hex(nData)
        else:
            if nData==2:
                aColor= ['red','blue']
            else:
                aColor=['red']

    if aHatch_in is not None:
        aHatch = aHatch_in
    else:
        aHatch= np.fill(nData, '+')

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )
    x = np.arange(  len(aLabel_x_in)  )
    dMin_x = -0.5
    dMax_x = len(aLabel_x_in)-0.5

     
   
    ax.set_ylabel(sLabel_y,fontsize=14)
    ax.set_title(sTitle,fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(aLabel_x_in)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)

    total_width = 0.6
    width = total_width / (nData-len(aReference_in))
    j = 0
    for i in np.arange(0, nData, 1):
        if i in aReference_in:
            pass 
        else:
            data1 = aData_in[i]
            
            rects = ax.bar( x - total_width * 0.5 + (j+0.5) * width, data1, width, label= aLabel_y_in[i], linestyle = aLinestyle_in[i],\
                                color = aColor[i], hatch = aHatch[i], edgecolor = "k")
            j=j+1
            pass

    for i in aReference_in:

        x0 = [-1, nData-len(aReference_in)]
        y0 = [aData_in[i][0], aData_in[i][0]]
        ax.plot( x0, y0, \
                 color = aColor[i], linestyle = 'dashed' ,\
                 marker = aMarker_in[i] ,\
                 label = aLabel_y_in[i])  

    if sLabel_info_in is not None:
        ax.text(0.1,0.9, sLabel_info_in, \
                        verticalalignment='center', horizontalalignment='left',\
                        transform=ax.transAxes, \
                        color='black', fontsize=13)

    if (iFlag_format_y ==1):
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter( sFormat_y ) )

    ax.set_xlim( dMin_x, dMax_x )
    ax.set_ylim( dMin_y, dMax_y )
    ax.grid(linewidth=1, color='gray', alpha=0.3, linestyle='--')
    ax.legend(bbox_to_anchor=aLocation_legend, \
              loc=sLocation_legend, \
              fontsize=15, \
              ncol= ncolumn)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()


    return
