import os, sys

import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

from pyearth.system.define_global_variables import *

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex



def barplot_data(aData_in,\
                 aLabel_x_in,\
                 aLabel_y_in,\
                 sFilename_out,\
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
                 sLabel_y_in = None, \
                 aLabel_legend_in = None,\
                 aLocation_legend_in =None,\
                 sFormat_y_in =None,\
                 sLocation_legend_in=None,\
                 sTitle_in = None):
    """
    Draw a barplot
    """

    aData_in = np.array(aData_in)
    pShape = aData_in.shape

    nData = pShape[0] #len(aLabel_x_in) this maybe checked 

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
        dMin_y = np.nanmin(aData_in) 

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
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(sLabel_y)
    ax.set_title(sTitle)
    ax.set_xticks(x)
    ax.set_xticklabels(aLabel_x_in)

    total_width = 0.7
    width = total_width / nData

    for i in np.arange(1, nData + 1, 1):

        data1 = aData_in[i-1]
        #need bettern control here, hatch (texture) can be used to betten represet some data
        if i == 1 or i ==2:
            rects = ax.bar( x - total_width * 0.5 + i * width, data1, width, label= aLabel_y_in[i-1],\
                            color = aColor[i-1], hatch = aHatch[i-1], edgecolor = "k")
        else:
            rects = ax.bar( x - total_width * 0.5 + i * width, data1, width, label= aLabel_y_in[i-1],\
                            color = aColor[i-1])#, hatch = aHatch[i-1], edgecolor = "k")
      
        pass

    if (iFlag_format_y ==1):
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter( sFormat_y ) )

    ax.set_ylim( dMin_y, dMax_y )
    ax.legend(bbox_to_anchor=aLocation_legend, \
              loc=sLocation_legend, \
              fontsize=12, \
              ncol= ncolumn)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()


    return
