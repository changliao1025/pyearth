import os, sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime


from pyearth.system.define_global_variables import *



def plot_time_series_data_multiple_temporal_resolution_fill(aTime_all, aData_all, \
                                  sFilename_out,\
                                  iDPI_in = None,\
                                  iFlag_trend_in = None, \
                                  iReverse_Y_in = None, \
                                  iSize_X_in = None, \
                                  iSize_Y_in = None, \
                                  dMax_Y_in =None, \
                                  dMin_Y_in = None, \
                                  dSpace_y_in=None,\
                                  aMarker_in =None,\
                                aColor_in =None,\
                                aLinestyle_in =None,\
                                  sLabel_Y_in = None, \
                                  aLabel_legend_in = None,\
                                  sTitle_in = None):

    nData = len(aData_all)
    aData = [aData_all]
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_trend_in is not None:
        iFlag_trend = 1
    else:
        iFlag_trend = 0

    if iReverse_Y_in is not None:
        iReverse_Y = 1
    else:
        iReverse_Y = 0

    if iSize_X_in is not None:
        iSize_X = iSize_X_in
    else:
        iSize_X = 12
    if iSize_Y_in is not None:
        iSize_Y = iSize_Y_in
    else:
        iSize_Y = 9

    if sLabel_Y_in is not None:
        sLabel_Y = sLabel_Y_in
    else:
        sLabel_Y = ''
    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = {}
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
        aColor= np.full(nData, 'black')
    if aLinestyle_in is not None:
        aLinestyle = aLinestyle_in
    else:

        aLinestyle = np.full(nData, '-')

    

    if dMax_Y_in is not None:
        dMax_Y = dMax_Y_in
    else:
        dMax_Y = np.nanmax(aData) * 1.2
    if dMin_Y_in is not None:
        dMin_Y = dMin_Y_in
    else:
        dMin_Y = np.nanmin(aData) #if it has negative value, change here
    if (dMax_Y <= dMin_Y ):
        return

    if dSpace_y_in is not None:        
        dSpace_y = dSpace_y_in
    else:       
        pass

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_X)
    fig.set_figheight( iSize_Y)
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )
    pYear = mdates.YearLocator(1)   # every year
    pMonth = mdates.MonthLocator()  # every month
    sYear_format = mdates.DateFormatter('%Y')
    
    for i in np.arange(1, nData+1):
        x1 = aTime_all[i-1]

        if i==1:
            #we need to plot low and high fill 
            y_top = aData_all[i-1][2]
            y_bot = aData_all[i-1][1]
            ax.fill_between(x1, y_top, y_bot,  facecolor='cornflowerblue')
            #plot mean
            y1 = aData_all[i-1][0]

            ax.plot( x1, y1, \
             color = aColor[i-1], linestyle = aLinestyle[i-1] ,\
             marker = aMarker[i-1] ,\
             label = aLabel_legend[i-1])

        else:
            y1 = aData_all[i-1]
        
            ax.plot( x1, y1, \
                 color = aColor[i-1], linestyle = aLinestyle[i-1] ,\
                 marker = aMarker[i-1] ,\
                 label = aLabel_legend[i-1])
    
        #calculate linear regression
        nan_index = np.where(y1 == missing_value)
        y1[nan_index] = np.nan
        good_index = np.where(  ~np.isnan(y1))
        if iFlag_trend ==1:
            x_dummy = np.array( [i.timestamp() for i in x1 ] )
            x_dummy = x_dummy[good_index]
            y_dummy = y1[good_index]
            coef = np.polyfit(x_dummy,y_dummy,1)
            poly1d_fn = np.poly1d(coef)
            mn = np.min(x_dummy)
            mx = np.max(x_dummy)
            x2 = [mn,mx]
            y2 = poly1d_fn(x2)
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
    ax.set_ylabel(sLabel_Y,fontsize=12)
    ax.set_title( sTitle, loc='center', fontsize=15)
    # round to nearest years...
    x_min = np.datetime64(aTime_all[1][0], 'Y') 
    x_max = np.datetime64(aTime_all[1][ -1 ], 'Y') + np.timedelta64(1, 'Y')
    ax.set_xlim(x_min, x_max)
    if dMax_Y < 1000 and dMax_Y > 0.001:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    else:
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
        
    
    ax.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y))
    
    if (iReverse_Y ==1):
        ax.set_ylim( dMax_Y, dMin_Y )
    else:
        ax.set_ylim( dMin_Y, dMax_Y )
    ax.legend(bbox_to_anchor=(1.0,1.0), loc="upper right", fontsize=12)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    #print('finished plotting')
