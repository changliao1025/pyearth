import os, sys
from datetime import datetime
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker


from pyearth.system.define_global_variables import *

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.choose_n_color import polylinear_gradient, rand_hex_color

def plot_time_series_data(aTime_all, \
                        aData_all, \
                          sFilename_out,\
                          iDPI_in = None,\
                          iFlag_log_in = None,\
                          iFlag_scientific_notation_in = None,\
                          ncolumn_in = None,\
                          aFlag_trend_in = None, \
                          iReverse_y_in = None, \
                          iSize_x_in = None, \
                          iSize_y_in = None, \
                          dMax_x_in = None, \
                          dMin_x_in = None, \
                          dMax_y_in = None, \
                          dMin_y_in = None, \
                          dSpace_y_in = None,\
                          aMarker_in = None,\
                          aColor_in = None,\
                          aLinestyle_in = None,\
                          sLabel_y_in = None, \
                          aLabel_legend_in = None,\
                          aLocation_legend_in =None,\
                          sDate_type_in = None,\
                          sFormat_y_in =None,\
                          sLocation_legend_in=None,\
                          sTitle_in = None):
    #find how many data will be plotted
    pShape = aData_all.shape
  
    nData = pShape[0]

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_log_in is not None:
        iFlag_log = 1
    else:
        iFlag_log = 0

    if iFlag_scientific_notation_in is not None:
        iFlag_scientific_notation = 1
    else:
        iFlag_scientific_notation = 0


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
            if (nData<=12):
                aColor= create_diverge_rgb_color_hex(nData)
            else:            
                a = rand_hex_color(num=2)
                b = polylinear_gradient(a, nData)
                aColor = b['hex']
                pass
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
        dMax_x = np.datetime64(np.nanmax(aTime_all), 'Y')
    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = np.datetime64(np.nanmin(aTime_all), 'Y')

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = np.nanmax(aData_all) #* 1.2

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aData_all) #if it has negative value, change here

    if (dMax_y <= dMin_y ):
        return

    if dSpace_y_in is not None:
        iFlag_space_y =1
        dSpace_y = dSpace_y_in
    else:
        iFlag_space_y =1
        dSpace_y = (dMax_y - dMin_y) /4.0
        pass

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )

    nYear = int( (dMax_x-dMin_x)/ 5 )
    pYear = mdates.YearLocator(nYear)   # every year
    pMonth = mdates.MonthLocator()  # every month

    if sDate_type_in is not None:
        if sDate_type_in == 'month':
            pMonth = mdates.MonthLocator(3)
        else:
            print(sDate_type_in)
    else:
        print(sDate_type_in)
        pass

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

    sYear_format = mdates.DateFormatter('%Y')

    #start loop for each data
    for i in np.arange(1, nData+1):

        x1 = aTime_all[i-1]
        y1 = aData_all[i-1]
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

    ax.xaxis.set_major_locator(pYear)
    ax.xaxis.set_minor_locator(pMonth)
    ax.xaxis.set_major_formatter(sYear_format)
    ax.tick_params(axis="x", labelsize=10)
    ax.tick_params(axis="y", labelsize=10)

    ax.set_xmargin(0.05)
    ax.set_ymargin(0.15)

    ax.set_xlabel('Year',fontsize=12)

    ax.set_title( sTitle, loc='center', fontsize=15)

    ax.set_xlim(dMin_x, dMax_x)

    #next Y axis
    ax.set_ylabel(sLabel_y,fontsize=12)

    if (iReverse_y ==1): #be careful here
        ax.set_ylim( dMax_y, dMin_y )
    else:
        ax.set_ylim( dMin_y, dMax_y )

    if iFlag_log ==1:
        #we need to change the ticklabel
        aLabel_y = []
        nlabel = int( (dMax_y- dMin_y) / dSpace_y) +1
        for i in np.arange( 0, nlabel, 1 ):
            ii = dMin_y + i * dSpace_y
            sTicklabel = r'$10^{{{}}}$'.format(ii)
            aLabel_y.append(sTicklabel)
            pass
        ticks = np.arange( 0, nlabel, 1 ) * dSpace_y + dMin_y
        ax.set_yticks( ticks)
        ax.set_yticklabels(aLabel_y)
        #if (iFlag_format_y ==1):
        #    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter( sFormat_y ) )

        
           
        #ax.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

        
        
    else:
        #not log

        if iFlag_scientific_notation ==1:
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1,1)) # you might need to change here
            ax.yaxis.set_major_formatter(formatter)
            #most time, when you use scientific notation, you may not need set the space,
            #but you may still set it using the method below

            pass
        else:
            if (iFlag_format_y ==1):
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter( sFormat_y ) )

            if (iFlag_space_y ==0):
                #ax.yaxis.set_major_locator(ticker.AutoLocator())
                #ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
                ax.yaxis.set_major_locator(ticker.MaxNLocator(prune='upper', nbins=5))
            else:
                ax.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y))
                ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

            pass

    if len(aLabel_legend[0]) > 0:
        ax.legend(bbox_to_anchor=aLocation_legend, \
              loc=sLocation_legend, \
              fontsize=12, \
              ncol= ncolumn)

    #save the result
    #plt.show()
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    #print('finished plotting')
