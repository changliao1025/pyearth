import os, sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime

from pyearth.system.define_global_variables import *



def plot_time_series_data_with_two_y_axis(aTime_all, aData_all, \
                                  sFilename_out,\
                                  iDPI_in = None,\
                                  iFlag_trend_in = None, \
                                  iReverse_y_in = None, \
                                  iSize_x_in = None, \
                                  iSize_y_in = None, \
                                  aMax_y_in =None, \
                                  aMin_y_in = None, \
                                  aSpace_y_in=None,\
                                  aMarker_in =None,\
                                aColor_in =None,\
                                aLinestyle_in =None,\
                                  aLabel_y_in = None, \
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

    if iReverse_y_in is not None:
        iReverse_y = 1
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

    if aLabel_y_in is not None:
        sLabel_y1 = aLabel_y_in[0]
        sLabel_y2 = aLabel_y_in[1]
    else:
        sLabel_y1 = ''
        sLabel_y2 = ''
    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
        sLabel_legend1 = aLabel_legend_in[0]
        sLabel_legend2 = aLabel_legend_in[1]
    else:
        sLabel_legend1 = ''
        sLabel_legend2 = ''
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
     

    if aMax_y_in is not None:
        dMax_y1 = aMax_y_in[0]
        dMax_y2 = aMax_y_in[1]
    else:
        dMax_y1  = np.nanmax(aData[0]) * 1.2
        dMax_y2  = np.nanmax(aData[1]) * 1.2
    if aMin_y_in is not None:
        dMin_y1 = aMin_y_in[0]
        dMin_y2 = aMin_y_in[1]
    else:
        #dMin_Y = np.nanmin(aData) #if it has negative value, change here
        dMin_y1  = np.nanmin(aData[0])
        dMin_y2  = np.nanmin(aData[1]) 

    #if (dMax_Y <= dMin_Y ):
    #    return

    if aSpace_y_in is not None:        
        dSpace_y1 = aSpace_y_in[0]
        dSpace_y2 = aSpace_y_in[1]
    else:       
        pass

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x)
    fig.set_figheight( iSize_y)
    ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4] )
    pYear = mdates.YearLocator(1)   # every year
    pMonth = mdates.MonthLocator()  # every month
    sYear_format = mdates.DateFormatter('%Y')
    #plot the first axis
    i=1
    x1 = aTime_all[i-1]
    
    y1 = aData_all[i-1]
    ax1.plot( x1, y1, \
            color = aColor[i-1], linestyle = aLinestyle[i-1] ,\
            marker = aMarker[i-1] , markersize =0.5,\
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
        ax1.plot(x2,y2, color = 'orange', linestyle = '-.',  linewidth=0.5)
    #plot the second axis
    ax2 = ax1.twinx()
    i=2
    x2 = aTime_all[i-1]
    y2 = aData_all[i-1]
    ax2.plot( x2, y2, \
            color = aColor[i-1], linestyle = aLinestyle[i-1] ,\
            marker = aMarker[i-1] ,\
            label = aLabel_legend[i-1])


    ax1.axis('on')
    ax1.grid(which='major', color='grey', linestyle='--', axis='y')
    
    #ax1.set_aspect(dRatio)  #this one set the y / x ratio
    ax1.xaxis.set_major_locator(pYear)
    ax1.xaxis.set_minor_locator(pMonth)
    ax1.xaxis.set_major_formatter(sYear_format)
    ax1.tick_params(axis="x", labelsize=10)
    ax1.tick_params(axis="y", labelsize=10)

    ax1.set_xmargin(0.05)
    ax1.set_ymargin(0.15)

    ax1.set_xlabel('Year',fontsize=12)
    ax1.set_ylabel(sLabel_y1,fontsize=12)
    ax2.set_ylabel(sLabel_y2,fontsize=12)
    ax1.set_title( sTitle, loc='center', fontsize=15)
    # round to nearest years...
    x_min = np.datetime64(aTime_all[1][0], 'Y') 
    x_max = np.datetime64(aTime_all[1][ -1 ], 'Y') + np.timedelta64(1, 'Y')
    ax1.set_xlim(x_min, x_max)
    if dMax_y1 < 1000 and dMax_y1 > 0.001:
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    else:
        ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
        
    
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y1))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(dSpace_y2))
    if (iReverse_y == 1):
        ax1.set_ylim( dMax_y1, dMin_y1 )
    else:
        ax1.set_ylim( dMin_y1, dMax_y1 )
        
    ax1.legend(bbox_to_anchor=(1.0,1.0), loc="upper right", fontsize=12)
    ax2.legend(bbox_to_anchor=(1.0,1.0), loc="upper right", fontsize=12)
    print(sFilename_out)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    #print('finished plotting')
