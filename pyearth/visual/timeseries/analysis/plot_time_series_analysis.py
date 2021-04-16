import os, sys
from datetime import datetime
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

import pandas as pd
from statsmodels.tsa.seasonal import STL
from statsmodels.tsa.stattools import adfuller

from pyearth.system.define_global_variables import *

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def plot_time_series_analysis(aTime, \
                              aData, \
                              sFilename_out,\
                              sVariable,\
                              iDPI_in = None,\
                              iFlag_without_raw_in= None,\
                              iFlag_log_in = None,\
                              iReverse_y_in = None, \
                              iSize_x_in = None, \
                              iSize_y_in = None, \
                              dMax_x_in = None, \
                              dMin_x_in = None, \
                              dMax_y_in = None, \
                              dMin_y_in = None, \
                              dSpace_x_in = None,\
                              dSpace_y_in = None,\
                              aMarker_in = None,\
                              aColor_in = None,\
                              aLinestyle_in = None,\
                              sLabel_x_in = None, \
                              sLabel_y_in = None, \
                              aLabel_legend_in = None,\
                              sDate_type_in = None,\
                              sFormat_y_in =None,\
                              sTitle_in = None):
    #find how many data will be plotted


    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_without_raw_in is not None:
        iFlag_without_raw = iFlag_without_raw_in
    else:
        iFlag_without_raw = 0

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
        sLabel_y = sVariable

    nData = 4
    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = np.array([ sVariable, 'Trend','Season','Residual' ])

    if sTitle_in is not None:
        sTitle = sTitle_in.title()
    else:
        sTitle = ''

    if aMarker_in is not None:
        aMarker = aMarker_in
    else:
        aMarker = ['o','.','*','+']

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
        aLinestyle =  ['-','--','-.' ,'solid']

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = np.datetime64(np.nanmax(aTime), 'Y')

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = np.datetime64(np.nanmin(aTime), 'Y')

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = np.nanmax(aData) * 1.2

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aData) * 0.8 #if it has negative value, change here
    if (dMax_y <= dMin_y ):
        return

    if dSpace_y_in is not None:
        iFlag_space_y =1
        dSpace_y = dSpace_y_in
    else:
        iFlag_space_y=0
        pass

    adf_test = adfuller(aData)
    #print(adf_test)
    print("ADF = " + str(adf_test[0]))
    print("p-value = " +str(adf_test[1])  )

    #fig = plt.figure( dpi=iDPI )
    #fig.set_figwidth( iSize_x )
    #fig.set_figheight( iSize_y )

    #pAxGrid = AxesGrid(fig, 111,
    #                    nrows_ncols=(4,1),
    #                    axes_pad=0.6,
    #                    label_mode='')  # note the empty label_mode


    pYear = mdates.YearLocator(1)   # every year
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

    sYear_format = mdates.DateFormatter('%Y')

    aData_tsa = pd.Series(aData, index=pd.date_range(aTime[0], \
                                                     periods=len(aTime), freq='M'), name = sVariable)

    
    #https://www.statsmodels.org/stable/generated/statsmodels.tsa.seasonal.STL.html#statsmodels.tsa.seasonal.STL
    stl = STL(aData_tsa, seasonal=13)
    aTSA = stl.fit()
    #part 1
    #plot time series



    if iFlag_without_raw == 1:
        aData_all = [aTSA.trend, aTSA.seasonal, aTSA.resid ]
        iSize_y = iSize_y * 0.75
        fig, pAxGrid = plt.subplots(nrows= len(aData_all), \
                                    figsize=(iSize_x, iSize_y),\
                                    sharex=True,  dpi=iDPI)
    else:
        aData_all = [aData, aTSA.trend, aTSA.seasonal, aTSA.resid ]
        fig, pAxGrid = plt.subplots(nrows= len(aData_all), \
                                    figsize=(iSize_x, iSize_y),\
                                    sharex=True,  dpi=iDPI)

    for i, ax in enumerate(pAxGrid):


        #ax.set_facecolor('#eafff5')
        if iFlag_without_raw == 1:
            ax.plot( aTime, aData_all[i], \
                     color = aColor[i+1], linestyle = aLinestyle[i+1] ,\
                     marker = aMarker[i+1] ,\
                     label = aLabel_legend[i+1],\
                     zorder=3)

        else:
            ax.plot( aTime, aData_all[i], \
                     color = aColor[i], linestyle = aLinestyle[i] ,\
                     marker = aMarker[i] ,\
                     label = aLabel_legend[i],\
                     zorder=3)

        if  iFlag_without_raw == 0:
            if i==0:
                if iReverse_y ==1:
                    ax.set_ylim( dMin_y,dMax_y  )

                ax.set_title(sTitle,fontsize=13)

                if iFlag_log_in ==1:
                    #we need to change the ticklabel
                    aLabel_y = []
                    for j in np.arange( dMin_y, dMax_y + 1, 1 ):
                        sTicklabel = r'$10^{{{}}}$'.format(int(j))
                        aLabel_y.append(sTicklabel)
                        pass

                    ax.set_yticks(np.arange( dMin_y, dMax_y +1, 1 ))
                    ax.set_yticklabels(aLabel_y)
                    pass
                pass
        else:
            if i==0:
                ax.set_title(sTitle,fontsize=13)


        if iFlag_without_raw == 1:
            ax.set_ylabel(aLabel_legend[i+1],fontsize=12)
            if i == 2:
                ax.plot((dMin_x, dMax_x), (0, 0), color='#000000', linestyle=':', zorder=2)

                ax.set_xlabel('Year',fontsize=12)

        else:
            ax.set_ylabel(aLabel_legend[i],fontsize=12)
            if i == 3:
                ax.plot((dMin_x, dMax_x), (0, 0), color='#000000', linestyle=':', zorder=2)
                ax.set_xlabel('Year',fontsize=12)


        ax.grid(which='major', color='lightgrey', linestyle=':', axis='y', zorder=1)


        ax.set_xlim(dMin_x, dMax_x)
        #ax.xaxis.set_major_locator(pYear)
        #ax.xaxis.set_minor_locator(pMonth)


    #we can now rewrite the plot function here





    #save the result
    #plt.show()
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    #print('finished plotting')
