
from datetime import datetime
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.choose_n_color import polylinear_gradient, rand_hex_color


def plot_time_series_data_w_variation(aTime_all,
                                      aData_all,
                                      aData_upper_all,
                                      aData_lower_all,
                                      sFilename_out,
                                      iDPI_in=None,
                                      iFlag_log_in=None,
                                      iFlag_miniplot_in=None,
                                      iFlag_scientific_notation_in=None,
                                      ncolumn_in=None,
                                      aFlag_trend_in=None,
                                      iReverse_y_in=None,
                                      iSize_x_in=None,
                                      iSize_y_in=None,
                                      dMax_x_in=None,
                                      dMin_x_in=None,
                                      dMax_y_in=None,
                                      dMin_y_in=None,
                                      dSpace_y_in=None,
                                      aMarker_in=None,
                                      aColor_in=None,
                                      aLinestyle_in=None,
                                      sLabel_y_in=None,
                                      aLabel_legend_in=None,
                                      aLocation_legend_in=None,
                                      sDate_type_in=None,
                                      sFormat_y_in=None,
                                      sLocation_legend_in=None,
                                      sTitle_in=None):
    # find how many data will be plotted
    aTime_all = np.array(aTime_all)
    aData_all = np.array(aData_all)
    pShape = aData_all.shape

    nData = pShape[0]

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_log_in is not None:
        iFlag_log = iFlag_log_in
    else:
        iFlag_log = 0

    if iFlag_scientific_notation_in is not None:
        iFlag_scientific_notation = iFlag_scientific_notation_in
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
        aLabel_legend = np.full(nData, '')

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    if aMarker_in is not None:
        aMarker = aMarker_in
    else:
        aMarker = np.full(nData, '+')

    if aColor_in is not None:
        aColor = aColor_in
        aColor_fill = aColor
    else:
        if (nData >= 3):
            if (nData <= 12):
                aColor = create_diverge_rgb_color_hex(nData)
                aColor_fill = aColor
            else:
                a = rand_hex_color(num=2)
                b = polylinear_gradient(a, nData)
                aColor = b['hex']
                aColor_fill = aColor
                pass
        else:
            if nData == 2:
                aColor = ['red', 'blue']
                aColor_fill = aColor
            else:
                aColor = ['navy']
                aColor_fill = ['lightblue']

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
        dMax_y = np.nanmax([aData_all, aData_upper_all, aData_lower_all])

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        # if it has negative value, change here
        dMin_y = np.nanmin([aData_all, aData_upper_all, aData_lower_all])

    if (dMax_y <= dMin_y):
        return
    else:
        iFlag_force_limit = 1
        if iFlag_force_limit == 1:
            pass
        else:
            dMin_y = dMin_y - 0.10 * (dMax_y-dMin_y)
            dMax_y = dMax_y + 0.25 * (dMax_y-dMin_y)  # leave space for legend

    if dSpace_y_in is not None:
        iFlag_space_y = 1
        dSpace_y = dSpace_y_in
    else:
        iFlag_space_y = 1
        dSpace_y = (dMax_y - dMin_y) / 4.0
        if dSpace_y < 1:
            pass
        else:
            dSpace_y = int(dSpace_y)
        pass

    if iFlag_miniplot_in is not None:
        iFlag_miniplot = iFlag_miniplot_in
        # set up location and range
        # dMin_mini_x = dMin_x + (dMax_x-dMin_x) * 0.7
        # dMax_mini_x =  dMax_x - (dMax_x-dMin_x) * 0.1
        # dMin_mini_y = dMin_y + (dMax_y-dMin_y) * 0.34
        # dMax_mini_y = dMax_y - (dMax_y-dMin_y) * 0.51

        dMin_mini_x = dMin_x + (dMax_x-dMin_x) * 0.7
        dMax_mini_x = dMax_x - (dMax_x-dMin_x) * 0.1
        dMin_mini_y = dMin_y + (dMax_y-dMin_y) * 0.05
        dMax_mini_y = dMax_y - (dMax_y-dMin_y) * 0.65
    else:
        iFlag_miniplot = 0

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.1, 0.8
    bottom, height = 0.1, 0.5
    dY_mini = 0.60
    dX_mini = 0.15
    width_mini = 0.28
    heigh_mini = 0.40
    rect_full = [left, bottom, width, height]
    rect_mini = [dY_mini, dX_mini, width_mini, heigh_mini]

    ax_full = plt.axes(rect_full)
    if iFlag_miniplot == 1:
        ax_mini = plt.axes(rect_mini)
        ax_all = [ax_full, ax_mini]
    else:
        ax_all = [ax_full]

    nYear = int((dMax_x-dMin_x) / 5)
    pYear = mpl.dates.YearLocator(nYear)   # every year
    pYear_min = mpl.dates.YearLocator(1)   # every year
    pMonth = mpl.dates.MonthLocator()  # every month
    pMonth_min = mpl.dates.MonthLocator(6)  # every month

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

    if ncolumn_in is not None:
        ncolumn = ncolumn_in
    else:
        ncolumn = 1

    sYear_format = mpl.dates.DateFormatter('%Y')

    # start loop for each data

    for iax in range(len(ax_all)):
        ax = ax_all[iax]
        ax.tick_params(direction='in', top=True, right=True)
        if iax == 0:
            aLegend_artist = []
            aLabel = []
        else:
            pass
        for i in np.arange(1, nData+1):
            x1 = aTime_all[i-1]
            y1 = aData_all[i-1]
            y1_upper = aData_upper_all[i-1]
            y1_lower = aData_lower_all[i-1]
            sc_var = ax.fill_between(x1,
                                     y1_upper,
                                     y1_lower,
                                     alpha=0.5, color=aColor_fill[i-1])

            tsp, = ax.plot(x1, y1,
                           color=aColor[i-1], linestyle=aLinestyle[i-1],
                           label=aLabel_legend[i-1])

            if iax == 0:
                aLegend_artist.append(tsp)
                aLabel.append(aLabel_legend[i-1])

                # calculate linear regression
                iFlag_trend = aFlag_trend[i-1]
                if iFlag_trend == 1:
                    nan_index = np.where(y1 == missing_value)
                    y1[nan_index] = np.nan
                    good_index = np.where(~np.isnan(y1))
                    x_dummy = np.array([i.timestamp() for i in x1])
                    x_dummy = x_dummy[good_index]
                    y_dummy = y1[good_index]
                    coef = np.polyfit(x_dummy, y_dummy, 1)
                    poly1d_fn = np.poly1d(coef)
                    mn = np.min(x_dummy)
                    mx = np.max(x_dummy)
                    x2 = [mn, mx]
                    y2 = poly1d_fn(x2)
                    x2 = [datetime.fromtimestamp(i) for i in x2]
                    ax.plot(x2, y2, color='orange',
                            linestyle='-.',  linewidth=0.5)

        # unqiue setting
        if iax == 0:
            ax.axis('on')
            ax.grid(which='major', color='grey', linestyle='--', axis='y')
            ax.xaxis.set_major_locator(pYear)
            ax.xaxis.set_minor_locator(pMonth)
            ax.xaxis.set_major_formatter(sYear_format)
            ax.tick_params(axis="x", labelsize=10)
            ax.tick_params(axis="y", labelsize=10)
            ax.set_xmargin(0.05)
            ax.set_ymargin(0.15)
            ax.set_xlabel('Year', fontsize=12)
            ax.set_title(sTitle, loc='center', fontsize=15)
            ax.set_xlim(dMin_x, dMax_x)
            ax.set_ylabel(sLabel_y, fontsize=12)

            if iFlag_log == 1:
                aLabel_y = list()
                if dSpace_y >= 1:
                    dSpace_y = int(dSpace_y)
                    nlabel = int((dMax_y - dMin_y) / dSpace_y) + 1
                    for i in np.arange(0, nlabel, 1):
                        ii = int(dMin_y) + i * dSpace_y
                        sTicklabel = r'$10^{{{}}}$'.format(int(ii))
                        aLabel_y.append(sTicklabel)
                        pass
                    ticks = np.arange(0, nlabel, 1) * dSpace_y + int(dMin_y)
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(aLabel_y)
                    pass
                else:
                    nlabel = int((dMax_y - dMin_y) / dSpace_y) + 1
                    for i in np.arange(0, nlabel, 1):
                        ii = int(dMin_y) + i * dSpace_y
                        iii = sFormat_y.format(ii)
                        sTicklabel = r'$10^{{{}}}$'.format(iii)
                        aLabel_y.append(sTicklabel)
                        pass
                    ticks = np.arange(0, nlabel, 1) * dSpace_y + dMin_y
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(aLabel_y)
                    pass

                ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
                pass
            else:
                if iFlag_scientific_notation == 1:
                    formatter = mpl.ticker.ScalarFormatter(useMathText=True)
                    formatter.set_scientific(True)
                    ax.yaxis.set_major_formatter(formatter)
                    pass
                else:
                    if (iFlag_space_y == 0):
                        ax.yaxis.set_major_locator(
                            mpl.ticker.MaxNLocator(prune='upper',    nbins=5))
                    else:
                        ax.yaxis.set_major_locator(
                            mpl.ticker.MultipleLocator(dSpace_y))
                        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
                        pass

                    if (iFlag_format_y == 1):
                        sFormat_y_dummy = sFormat_y.replace("{", "{x")
                        ax.yaxis.set_major_formatter(
                            mpl.ticker.StrMethodFormatter(sFormat_y_dummy))
                        pass
                    else:
                        pass

                    pass

            if (iReverse_y == 1):  # be careful here
                ax.set_ylim(dMax_y, dMin_y)
            else:
                ax.set_ylim(dMin_y, dMax_y)

            ax.legend(aLegend_artist, aLabel, bbox_to_anchor=aLocation_legend,
                      loc=sLocation_legend, fontsize=12, ncol=ncolumn)
        else:
            if iFlag_log == 1:
                aLabel_y = list()
                dSpace_y_mini = dSpace_y / 3
                if dSpace_y_mini >= 1:
                    dSpace_y_mini = int(dSpace_y_mini)
                    nlabel = int((dMax_y - dMin_y) / dSpace_y_mini) + 1
                    for i in np.arange(0, nlabel, 1):
                        ii = int(dMin_y) + i * dSpace_y_mini
                        sTicklabel = r'$10^{{{}}}$'.format(int(ii))
                        aLabel_y.append(sTicklabel)
                        pass
                    ticks = np.arange(0, nlabel, 1) * \
                        dSpace_y_mini + int(dMin_y)
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(aLabel_y)
                else:
                    nlabel = int((dMax_y - dMin_y) / dSpace_y_mini) + 1
                    for i in np.arange(0, nlabel, 1):
                        ii = int(dMin_y) + i * dSpace_y_mini
                        iii = sFormat_y.format(ii)
                        sTicklabel = r'$10^{{{}}}$'.format(iii)
                        aLabel_y.append(sTicklabel)
                        pass
                    ticks = np.arange(0, nlabel, 1) * dSpace_y_mini + dMin_y
                    ax.set_yticks(ticks)
                    ax.set_yticklabels(aLabel_y)
                    pass
                ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
                pass
            else:
                if iFlag_scientific_notation == 1:
                    formatter = mpl.ticker.ScalarFormatter(useMathText=True)
                    formatter.set_scientific(True)
                    ax.yaxis.set_major_formatter(formatter)
                    pass
                else:
                    dSpace_y_mini = dSpace_y / 4
                    ax.yaxis.set_major_locator(
                        mpl.ticker.MultipleLocator(dSpace_y_mini))
                    ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
                    if (iFlag_format_y == 1):
                        sFormat_y_dummy = sFormat_y.replace("{", "{x")
                        ax.yaxis.set_major_formatter(
                            mpl.ticker.StrMethodFormatter(sFormat_y_dummy))
                    pass
            ax.xaxis.set_major_locator(pYear_min)
            ax.xaxis.set_minor_locator(pMonth_min)
            ax.xaxis.set_major_formatter(sYear_format)
            ax.set_xlim(dMin_mini_x, dMax_mini_x)
            if (iReverse_y == 1):
                ax.set_ylim(dMax_mini_y, dMin_mini_y)
            else:
                ax.set_ylim(dMin_mini_y, dMax_mini_y)
            pass

    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close('all')
    plt.clf()
