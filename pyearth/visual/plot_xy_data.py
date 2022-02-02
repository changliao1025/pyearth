import os, sys
from datetime import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def plot_xy_data(aX_all,  aY_all,  sFilename_out,iDPI_in = None, aFlag_trend_in = None,  iReverse_y_in = None,  iSize_x_in = None,  iSize_y_in = None,  ncolumn_in = None,  dMax_x_in = None,  dMin_x_in = None,  dMax_y_in = None, dMin_y_in = None,  dSpace_y_in = None, aMarker_in = None, aColor_in = None, aLinestyle_in = None,  aLabel_point_in = None,  aTick_label_x_in = None,  aLocation_legend_in =None, sLabel_x_in = None, sLabel_y_in = None,  aLabel_legend_in = None, sLocation_legend_in=None, sFormat_y_in =None, sTitle_in = None):
    """
    Plot a X-Y data

    Args:
        aX_all ([type]): [description]
        aY_all ([type]): [description]
        sFilename_out ([type]): [description]
        iDPI_in ([type], optional): [description]. Defaults to None.
        aFlag_trend_in ([type], optional): [description]. Defaults to None.
        iReverse_y_in ([type], optional): [description]. Defaults to None.
        iSize_x_in ([type], optional): [description]. Defaults to None.
        iSize_y_in ([type], optional): [description]. Defaults to None.
        ncolumn_in ([type], optional): [description]. Defaults to None.
        dMax_x_in ([type], optional): [description]. Defaults to None.
        dMin_x_in ([type], optional): [description]. Defaults to None.
        dMax_y_in ([type], optional): [description]. Defaults to None.
        dMin_y_in ([type], optional): [description]. Defaults to None.
        dSpace_y_in ([type], optional): [description]. Defaults to None.
        aMarker_in ([type], optional): [description]. Defaults to None.
        aColor_in ([type], optional): [description]. Defaults to None.
        aLinestyle_in ([type], optional): [description]. Defaults to None.
        aLabel_point_in ([type], optional): [description]. Defaults to None.
        aTick_label_x_in ([type], optional): [description]. Defaults to None.
        aLocation_legend_in ([type], optional): [description]. Defaults to None.
        sLabel_x_in ([type], optional): [description]. Defaults to None.
        sLabel_y_in ([type], optional): [description]. Defaults to None.
        aLabel_legend_in ([type], optional): [description]. Defaults to None.
        sLocation_legend_in ([type], optional): [description]. Defaults to None.
        sFormat_y_in ([type], optional): [description]. Defaults to None.
        sTitle_in ([type], optional): [description]. Defaults to None.
    """

    #find how many data will be plotted
    nData = len(aY_all)

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if aFlag_trend_in is not None:
        aFlag_trend = aFlag_trend_in
    else:
        aFlag_trend = np.full(nData, 0)

    if iReverse_y_in is not None:
        iReverse_y = iReverse_y_in
    else:
        iReverse_y = 0

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 12

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 9
    
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

    if sLabel_x_in is not None:
        sLabel_x = sLabel_x_in
    else:
        sLabel_x = ''

    if sLabel_y_in is not None:
        sLabel_y = sLabel_y_in
    else:
        sLabel_y = ''

    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = np.full(nData,'')

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    if aMarker_in is not None:
        aMarker = aMarker_in
    else:
        aMarker=np.full(nData, '+')

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

    if aLinestyle_in is not None:
        aLinestyle = aLinestyle_in
    else:
        aLinestyle = np.full(nData, '-')

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = np.nanmax(aX_all)
    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = np.nanmin(aX_all)

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = np.nanmax(aY_all) * 1.0
    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aY_all) #if it has negative value, change here
    if (dMax_y <= dMin_y ):
        return

    if dSpace_y_in is not None:
        iFlag_space_y =1
        dSpace_y = dSpace_y_in
    else:
        iFlag_space_y=0
        pass

    if aLabel_point_in is not None:
        iFlag_label =1
        aLabel_point = aLabel_point_in
        pass
    else:
        iFlag_label=0
        pass

    if aTick_label_x_in is not None:
        iFlag_replace_xtick = 1
        xtick_labels = aTick_label_x_in
        pass
    else:
        iFlag_replace_xtick = 0
        pass


    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )

  
    if sFormat_y_in is not None:
        iFlag_format_y = 1
        sFormat_y = sFormat_y_in
    else:
        iFlag_format_y = 0
 

    #start loop for each data
    for i in np.arange(1, nData+1):

        x1 = aX_all[i-1]
        y1 = aY_all[i-1]
        ax.plot( x1, y1, \
                 color = aColor[i-1], linestyle = aLinestyle[i-1] ,\
                 marker = aMarker[i-1] ,\
                 label = aLabel_legend[i-1])

        #calculate linear regression
        iFlag_trend = aFlag_trend[i-1]
        if iFlag_trend ==1:
            nan_index = np.where(y1 == missing_value)
            y1[nan_index] = np.nan
            good_index = np.where(  ~np.isnan(y1))

            x_dummy = np.array( [i.timestamp() for i in x1 ] )
            x_dummy = x_dummy[good_index]
            y_dummy = y1[good_index]
            coef = np.polyfit(x_dummy,y_dummy,1)
            poly1d_fn = np.poly1d(coef)
            mn=np.min(x_dummy)
            mx=np.max(x_dummy)
            x2=[mn,mx]
            y2=poly1d_fn(x2)
            x2 = [datetime.fromtimestamp(i) for i in x2 ]
            ax.plot(x2,y2, color = 'orange', linestyle = '-.',  linewidth=0.5)

    ax.axis('on')
    ax.grid(which='major', color='grey', linestyle='--', axis='y')

    #ax.set_aspect(dRatio)  #this one set the y / x ratio

   
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)

    ax.set_xmargin(0.05)
    ax.set_ymargin(0.15)

    ax.set_xlabel(sLabel_x,fontsize=12)
    ax.set_ylabel(sLabel_y,fontsize=12)
    ax.set_title( sTitle, loc='center', fontsize=15)

    ax.set_xlim(dMin_x, dMax_x)

    #if dMax_y < 1000 and dMax_y > 0.001:
    #    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    #else:
    #    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
    if (iFlag_format_y ==1):
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(  sFormat_y ))

    if (iFlag_space_y ==0):
        ax.yaxis.set_major_locator(ticker.AutoLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    else:
        ax.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y))

    if (iReverse_y ==1):
        ax.set_ylim( dMax_y, dMin_y )
    else:
        ax.set_ylim( dMin_y, dMax_y )
        pass

    
    
        

    if iFlag_label ==1:
        aText = []
        for i in np.arange(1, len(aLabel_point)+1, 1):

            
            sLabel_point= aLabel_point[i-1]
            aText.append(ax.text(aX_all[1][i-1], aY_all[1][i-1], sLabel_point, color = 'red'))
            
            pass

        adjust_text(aText, \
            only_move={'points':'y', 'texts':'xy'}, \
            arrowprops=dict(arrowstyle="->", color='r', lw=0.5) )

        pass

    if iFlag_replace_xtick ==1:
        ax.set_xticks(aX_all[0])
        ax.set_xticklabels(xtick_labels,fontsize= 8 )
        pass

    ax.legend(bbox_to_anchor=aLocation_legend, \
              loc=sLocation_legend, \
              fontsize=8, \
              ncol= ncolumn)

   
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
 
