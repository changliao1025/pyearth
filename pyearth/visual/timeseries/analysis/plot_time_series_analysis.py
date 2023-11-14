import os
import sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def plot_time_series_analysis(aTime,
                              aData,
                              sFilename_out,
                              sVariable,
                              iDPI_in=None,
                              iFlag_without_raw_in=None,
                              iFlag_log_in=None,
                              iFlag_scientific_notation_in=None,
                              iReverse_y_in=None,
                              iSize_x_in=None,
                              iSize_y_in=None,
                              dMax_x_in=None,
                              dMin_x_in=None,
                              dMax_y_in=None,
                              dMin_y_in=None,
                              dSpace_x_in=None,
                              dSpace_y_in=None,
                              aMarker_in=None,
                              aColor_in=None,
                              aLinestyle_in=None,
                              sLabel_x_in=None,
                              sLabel_y_in=None,
                              aLabel_legend_in=None,
                              aLocation_legend_in=None,
                              sLocation_legend_in=None,
                              sDate_type_in=None,
                              sFormat_y_in=None,
                              sTitle_in=None):

    try:
        import pandas as pd
        from statsmodels.tsa.seasonal import STL
        from statsmodels.tsa.stattools import adfuller
    except ImportError as e:
        raise ImportError(
            "The package 'pandas and statsmodels' is required for this function to run.") from e

    aTime = np.array(aTime)
    aData = np.array(aData)
    pShape = aData.shape

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
        sLabel_y = ''

    nData = 4
    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = np.full(nData, '')

    if sTitle_in is not None:
        sTitle = sTitle_in.capitalize()  # title will also convert units
    else:
        sTitle = ''

    if aMarker_in is not None:
        aMarker = aMarker_in
    else:
        aMarker = ['o', '.', '*', '+']

    if aColor_in is not None:
        aColor = aColor_in
    else:
        aColor = create_diverge_rgb_color_hex(nData)

    if aLinestyle_in is not None:
        aLinestyle = aLinestyle_in
    else:
        aLinestyle = ['-', '--', '-.', 'solid']

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
        dMax_y = np.nanmax(aData)

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aData)

    if (dMax_y <= dMin_y):
        return
    else:
        dMin_y = dMin_y - 0.13 * (dMax_y-dMin_y)
        dMax_y = dMax_y + 0.13 * (dMax_y-dMin_y)

    if dSpace_y_in is not None:
        iFlag_space_y = 1
        dSpace_y = dSpace_y_in
    else:
        iFlag_space_y = 0
        dSpace_y = (dMax_y - dMin_y) / 4.0
        if dSpace_y < 1:
            pass
        else:
            dSpace_y = int(dSpace_y)
        pass

    adf_test = adfuller(aData)
    print("ADF = " + str(adf_test[0]))
    print("p-value = " + str(adf_test[1]))
    print("Critical value = ", (adf_test[4]))

    pYear = mpl.dates.YearLocator(1)   # every year
    pMonth = mpl.dates.MonthLocator()  # every month
    if sDate_type_in is not None:
        if sDate_type_in == 'month':
            pMonth = mpl.dates.MonthLocator(3)
        else:
            print(sDate_type_in)
    else:
        pass

    if sFormat_y_in is not None:
        iFlag_format_y = 1
        sFormat_y = sFormat_y_in
    else:
        iFlag_format_y = 0
        sFormat_y = '{:.1f}'

    if sLocation_legend_in is not None:
        sLocation_legend = sLocation_legend_in
    else:
        sLocation_legend = "upper right"

    if aLocation_legend_in is not None:
        aLocation_legend = aLocation_legend_in
    else:
        aLocation_legend = (1.0, 1.0)

    sYear_format = mpl.dates.DateFormatter('%Y')

    aTS = pd.date_range(aTime[0],         periods=len(aTime),      freq='M')
    aData_tsa = pd.Series(aData,    index=aTS,      name=sVariable)

    # https://www.statsmodels.org/stable/generated/statsmodels.tsa.seasonal.STL.html#statsmodels.tsa.seasonal.STL
    stl = STL(aData_tsa, seasonal=13)
    aTSA = stl.fit()
    # part 1
    # plot time series
    if iFlag_without_raw == 1:
        aLabel_y = np.array(['Trend', 'Season', 'Residual'])
        aData_all = [aTSA.trend, aTSA.seasonal, aTSA.resid]
        iSize_y = iSize_y * 0.75
        fig, pAxGrid = plt.subplots(nrows=len(aData_all),
                                           figsize=(iSize_x, iSize_y),
                                           sharex=True,  dpi=iDPI)

        for i, ax in enumerate(pAxGrid):
            tsp, = ax.plot(aTime, aData_all[i],
                           color=aColor[i+1], linestyle=aLinestyle[i+1],
                           marker=aMarker[i+1],
                           zorder=3)
            if i == 0:
                ax.set_title(sTitle, fontsize=13)

                # we will not use legend method because it has 3 panels
                if aLabel_legend_in is not None:
                    # plot the first on the top
                    sText = aLabel_legend_in[0]
                    dLocation = 0.95
                    ax.text(0.03, dLocation, sText,
                            verticalalignment='top', horizontalalignment='left',
                            transform=ax.transAxes,
                            color='black', fontsize=10)
                    # plot the remaining on the bot
                    nlegend = len(aLabel_legend_in)
                    for i in range(1, nlegend, 1):
                        sText = aLabel_legend_in[i]
                        dLocation = nlegend * 0.06 - i * 0.05 - 0.03
                        ax.text(0.03, dLocation, sText,
                                verticalalignment='top', horizontalalignment='left',
                                transform=ax.transAxes,
                                color='black', fontsize=10)
                        pass

            if i == 2:
                ax.plot((dMin_x, dMax_x), (0, 0),
                        color='#000000', linestyle=':', zorder=2)
                ax.set_xlabel('Year', fontsize=12)
            ax.set_ylabel(aLabel_y[i], fontsize=12)
            ax.grid(which='major', color='lightgrey',
                    linestyle=':', axis='y', zorder=1)
            ax.set_xlim(dMin_x, dMax_x)
            ax.xaxis.set_major_locator(pYear)
            ax.xaxis.set_minor_locator(pMonth)
            ax.xaxis.set_major_formatter(sYear_format)

    else:
        # include raw data
        aLabel_y = np.array([sLabel_y, 'Trend', 'Season', 'Residual'])
        aData_all = [aData, aTSA.trend, aTSA.seasonal, aTSA.resid]
        fig, pAxGrid = plt.subplots(nrows=len(aData_all),
                                           figsize=(iSize_x, iSize_y),
                                           sharex=True,  dpi=iDPI)
        for i, ax in enumerate(pAxGrid):
            tsp, = ax.plot(aTime, aData_all[i],
                           color=aColor[i], linestyle=aLinestyle[i],
                           marker=aMarker[i],
                           zorder=3)
            # the first plot has title
            if i == 0:
                # might have to consider log label here, refer to time series plot
                aTickLabel_y = list()
                ax.set_title(sTitle, fontsize=13)
                if iReverse_y == 1:
                    ax.set_ylim(dMin_y, dMax_y)

            # the bottom plot has x label and provided y label
            if i == 3:
                ax.plot((dMin_x, dMax_x), (0, 0),
                        color='#000000', linestyle=':', zorder=2)
                ax.set_xlabel('Year', fontsize=12)

            ax.set_ylabel(aLabel_y[i], fontsize=12)
            ax.grid(which='major', color='lightgrey',
                    linestyle=':', axis='y', zorder=1)
            ax.set_xlim(dMin_x, dMax_x)
            ax.xaxis.set_major_locator(pYear)
            ax.xaxis.set_minor_locator(pMonth)
            ax.xaxis.set_major_formatter(sYear_format)

    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close('all')
    plt.clf()
