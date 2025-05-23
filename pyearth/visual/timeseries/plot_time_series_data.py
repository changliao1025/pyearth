
from datetime import datetime
import matplotlib.dates as mdates
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch

from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.choose_n_color import polylinear_gradient, rand_hex_color
from pyearth.visual.formatter import log_formatter, MathTextSciFormatter


def plot_time_series_data(aTime_all,
                          aData_all,
                          sFilename_out=None,
                          iDPI_in=None,
                          iFlag_log_in=None,
                          iFlag_scientific_notation_in=None,
                          iFlag_miniplot_in=None,
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
                          aLocation_miniplot_in = None,
                          aLocation_legend_in=None,
                          sDate_type_in=None,
                          sFormat_y_in=None,
                          sFont_in=None,
                          sLocation_legend_in=None,
                          sTitle_in=None):
    """
    Plot time series data

    Args:
        aTime_all (list): List of datatime for each date series, each data time sereis can be different
        aData_all (list): List of data for each date series
        sFilename_out (str): The output figure filename
        iDPI_in (_type_, optional): _description_. Defaults to None.
        iFlag_log_in (_type_, optional): _description_. Defaults to None.
        iFlag_scientific_notation_in (_type_, optional): _description_. Defaults to None.
        iFlag_miniplot_in (_type_, optional): _description_. Defaults to None.
        ncolumn_in (_type_, optional): _description_. Defaults to None.
        aFlag_trend_in (_type_, optional): _description_. Defaults to None.
        iReverse_y_in (_type_, optional): _description_. Defaults to None.
        iSize_x_in (_type_, optional): _description_. Defaults to None.
        iSize_y_in (_type_, optional): _description_. Defaults to None.
        dMax_x_in (_type_, optional): _description_. Defaults to None.
        dMin_x_in (_type_, optional): _description_. Defaults to None.
        dMax_y_in (_type_, optional): _description_. Defaults to None.
        dMin_y_in (_type_, optional): _description_. Defaults to None.
        dSpace_y_in (_type_, optional): _description_. Defaults to None.
        aMarker_in (_type_, optional): _description_. Defaults to None.
        aColor_in (_type_, optional): _description_. Defaults to None.
        aLinestyle_in (_type_, optional): _description_. Defaults to None.
        sLabel_y_in (_type_, optional): _description_. Defaults to None.
        aLabel_legend_in (_type_, optional): _description_. Defaults to None.
        aLocation_legend_in (_type_, optional): _description_. Defaults to None.
        sDate_type_in (_type_, optional): _description_. Defaults to None.
        sFormat_y_in (_type_, optional): _description_. Defaults to None.
        sLocation_legend_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
    """

    aTime_all = np.array(aTime_all)
    # each list is a data series, but length may be different
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

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    if aMarker_in is not None:
        aMarker = aMarker_in
    else:
        aMarker = np.full(nData, '+')

    if aColor_in is not None:
        aColor = aColor_in
    else:
        if (nData >= 3):
            if (nData <= 12):
                aColor = create_diverge_rgb_color_hex(nData)
            else:
                a = rand_hex_color(num=2)
                b = polylinear_gradient(a, nData)
                aColor = b['hex']
                pass
        else:
            if nData == 2:
                aColor = ['red', 'blue']
            else:
                aColor = ['red']

    if aLinestyle_in is not None:
        aLinestyle = aLinestyle_in
    else:
        aLinestyle = np.full(nData, '-')

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = datetime(1, 1, 1)
        for i in range(nData):
            dummy = np.nanmax(aTime_all[i])
            dMax_x = np.max([dMax_x, dummy])
        dMax_x_y = np.datetime64(dMax_x, 'Y')
        dMax_x_m = np.datetime64(dMax_x, 'M')

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        # dMin_x = np.datetime64(np.nanmin(aTime_all), 'Y')
        dMin_x = datetime(9999, 1, 1)
        for i in range(nData):
            dummy = np.nanmin(aTime_all[i])
            dMin_x = np.min([dMin_x, dummy])
        dMin_x_m = np.datetime64(dMin_x, 'M')
        dMin_x_y = np.datetime64(dMin_x, 'Y')

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
        iFlag_force_limit_y1 = 1
    else:
        iFlag_force_limit_y1 = 0
        dMax_y = np.nanmax(aData_all[0])
        for i in range(1, nData):
            dummy = np.nanmax(aData_all[i])
            dMax_y = np.max([dMax_y, dummy])


    if dMin_y_in is not None:
        dMin_y = dMin_y_in
        iFlag_force_limit_y0 = 1
    else:
        iFlag_force_limit_y0 = 0
        dMin_y = np.nanmin(aData_all[0])
        for i in range(1, nData):
            dummy = np.nanmin(aData_all[i])
            dMin_y = np.min([dMin_y, dummy])

    if (dMax_y <= dMin_y):
        return
    else:
        if iFlag_force_limit_y0 != 1:
            dMin_y = dMin_y - 0.10 * (dMax_y-dMin_y)
        if iFlag_force_limit_y1 != 1:
            dMax_y = dMax_y + 0.10 * (dMax_y-dMin_y)

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

    nYear = int(dMax_x_y-dMin_x_y) + 1
    if nYear > 3:
        dMin_x = dMin_x_y
        dMax_x = dMax_x_y
    else:
        dMin_x = dMin_x_m
        dMax_x = dMax_x_m

    if iFlag_miniplot_in is not None:
        iFlag_miniplot = iFlag_miniplot_in
        # set up location and range
        dMin_mini_x = dMin_x + (dMax_x-dMin_x) * 0.7
        dMax_mini_x = dMax_x - (dMax_x-dMin_x) * 0.1
        dMin_mini_y = dMin_y + (dMax_y-dMin_y) * 0.1
        dMax_mini_y = dMin_y + (dMax_y-dMin_y) * 0.4

        #add a rectangle to the mini plot
        # Convert datetime values to numerical format
        dMin_mini_x_num = mdates.date2num(dMin_mini_x)
        dMax_mini_x_num = mdates.date2num(dMax_mini_x)

        # Add a rectangle to the main plot
        rect = patches.Rectangle(
        (dMin_mini_x_num, dMin_mini_y),  # Bottom-left corner (x, y)
        dMax_mini_x_num - dMin_mini_x_num,  # Width
        dMax_mini_y - dMin_mini_y,  # Height
        linewidth=2,  # Border thickness
        edgecolor='black',  # Border color
        facecolor='none',  # Transparent fill
        linestyle='--'  # Dashed border
        )

    else:
        iFlag_miniplot = 0

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.1, 0.8

    if iFlag_miniplot == 1:
        bottom, height = 0.55, 0.45
        if aLocation_miniplot_in is not None:
            rect_mini = aLocation_miniplot_in
        else:
            dX_mini = left #0.60
            dY_mini = 0.0 #0.20
            width_mini = width #0.28
            heigh_mini = height # 0.30
            rect_mini = [dX_mini, dY_mini, width_mini, heigh_mini]
            pass

        rect_full = [left, bottom, width, height]
        ax_full = plt.axes(rect_full)
        ax_mini = plt.axes(rect_mini)
        ax_all = [ax_full, ax_mini]

    else:
        bottom, height = 0.1, 0.5
        rect_full = [left, bottom, width, height]
        ax_full = plt.axes(rect_full)
        ax_all = [ax_full]

    if nYear <= 3:
        pYear = mpl.dates.YearLocator()   # every year
        pMonth = mpl.dates.MonthLocator()  # every month

        pYear_min = mpl.dates.YearLocator(1)   # every year
        pMonth_min = mpl.dates.MonthLocator(3)  # every 3 month
    else:
        pYear = mpl.dates.YearLocator(2)   # every 2 year
        pMonth = mpl.dates.MonthLocator(6)  # every month

        pYear_min = mpl.dates.YearLocator(1)   # every year
        pMonth_min = mpl.dates.MonthLocator(6)  # every 6 month

    if sDate_type_in is not None:
        if sDate_type_in == 'month':
            pMonth = mpl.dates.MonthLocator(3)
        else:
            pass
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
        sLocation_legend = 'best'  # remove the setting so it becomes automatical

    if aLocation_legend_in is not None:
        aLocation_legend = aLocation_legend_in
    else:
        aLocation_legend = None  # (1.0,1.0)

    if ncolumn_in is not None:
        ncolumn = ncolumn_in
    else:
        ncolumn = 1

    sYear_format = mpl.dates.DateFormatter('%Y')
    sMonth_format = mpl.dates.DateFormatter('%Y-%m')
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
            tsp, = ax.plot(x1, y1,
                           color=aColor[i-1], linestyle=aLinestyle[i-1],
                           marker=aMarker[i-1])

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

            else:
                pass

        # unqiue setting
        if iax == 0:
            ax.axis('on')
            ax.grid(which='major', color='grey', linestyle='--', axis='y')
            if nYear <= 3:
                ax.xaxis.set_major_locator(pYear)
                ax.xaxis.set_minor_locator(pMonth)
                ax.xaxis.set_major_formatter(sMonth_format)
                #set tick for each month
                ax.xaxis.set_minor_formatter(mpl.dates.DateFormatter('%m'))
                #also set the label for tick
                # Convert tick positions to datetime before formatting
                ax.set_xticklabels([mpl.dates.num2date(x).strftime('%Y-%m') for x in ax.get_xticks()])
            else:
                ax.xaxis.set_major_locator(pYear)
                ax.xaxis.set_minor_locator(pMonth)
                ax.xaxis.set_major_formatter(sYear_format)
            ax.tick_params(axis="x", labelsize=10)
            ax.tick_params(axis="y", labelsize=10)
            ax.set_xmargin(0.05)
            ax.set_ymargin(0.15)
            if (iReverse_y == 1):
                ax.set_ylim(dMax_y, dMin_y)
            else:
                ax.set_ylim(dMin_y, dMax_y)
            # y axis labels are different because space is different

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
                    # formatter = mpl.ticker.ScalarFormatter(useMathText=True)
                    # formatter.set_scientific(True)
                    # y0 = int(np.log10(dMin_y))
                    # y1=  int(np.log10(dMax_y))
                    # formatter.set_powerlimits(( y0, y1))
                    # ax.yaxis.set_major_formatter(formatter)
                    ax.yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))

                    pass
                else:
                    if (iFlag_space_y == 0):
                        ax.yaxis.set_major_locator(
                            mpl.ticker.MaxNLocator(prune='upper', nbins=5))
                    else:
                        ax.yaxis.set_major_locator(
                            mpl.ticker.MultipleLocator(dSpace_y))
                        ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())
                    if (iFlag_format_y == 1):
                        sFormat_y_dummy = sFormat_y.replace("{", "{x")
                        ax.yaxis.set_major_formatter(
                            mpl.ticker.StrMethodFormatter(sFormat_y_dummy))
                    pass
            ax.set_title(sTitle, loc='center', fontsize=15)
            ax.set_xlim(dMin_x_m, dMax_x_m)
            ax.set_ylabel(sLabel_y, fontsize=12)

            ax.set_xlabel('Year', fontsize=12)
            ax.legend(aLegend_artist, aLabel, bbox_to_anchor=aLocation_legend,
                      loc=sLocation_legend, fontsize=10, ncol=ncolumn )
            pass
            if iFlag_miniplot == 1:
                rect_bottom_left = ax.transData.transform((dMin_mini_x_num, dMin_mini_y))
                rect_top_right = ax.transData.transform((dMax_mini_x_num, dMax_mini_y))
                # Convert display coordinates to figure coordinates (normalized [0, 1] space)
                rect_bottom_left_fig = fig.transFigure.inverted().transform(rect_bottom_left)
                rect_top_right_fig = fig.transFigure.inverted().transform(rect_top_right)
                # Extract normalized coordinates
                dMin_mini_x_fig, dMin_mini_y_fig = rect_bottom_left_fig
                dMax_mini_x_fig, dMax_mini_y_fig = rect_top_right_fig

                #add a label at upper left corner
                sLabel_info_in= '(a)'
                ax.text(0.05, 0.9, sLabel_info_in,
                    verticalalignment='center', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=13)

        else:
            if iFlag_log == 1:
                aLabel_y = list()
                dSpace_y_mini = dSpace_y / 4
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

            bbox_mini = ax.get_position()
            left_x_mini = bbox_mini.x0   # Center x of the miniplot
            right_x_mini = bbox_mini.x0 + bbox_mini.width   # Center x of the miniplot
            top_y_mini = bbox_mini.y0 + bbox_mini.height   # Center y of the miniplot
            sLabel_info_in= '(b)'
            ax.text(0.05, 0.9, sLabel_info_in,
                    verticalalignment='center', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black', fontsize=13)

    if iFlag_miniplot == 1:
        # Define the arrows pointing to the miniplot
        arrow_left = FancyArrowPatch(
            (dMin_mini_x_fig, dMin_mini_y_fig),  # Start point (bottom-left corner of the rectangle)
            (left_x_mini, top_y_mini),  # End point (center of the miniplot)
            arrowstyle="->",  # Arrow style
            color="black",  # Arrow color
            mutation_scale=15,  # Scale of the arrow
            linewidth=1.5,  # Arrow thickness
            clip_on=False,
            transform=fig.transFigure )
        arrow_right = FancyArrowPatch(
            (dMax_mini_x_fig, dMin_mini_y_fig),  # Start point (bottom-right corner of the rectangle)
            (right_x_mini, top_y_mini),  # End point (center of the miniplot)
            arrowstyle="->",  # Arrow style
            color="black",  # Arrow color
            mutation_scale=15,  # Scale of the arrow
            linewidth=1.5,  # Arrow thickness
            clip_on=False,
            transform=fig.transFigure  )
        ax = ax_all[0]
        ax.add_patch(rect)
        ax.add_patch(arrow_left)
        ax.add_patch(arrow_right)


    if sFilename_out is not None:
        plt.savefig(sFilename_out, bbox_inches='tight')
        plt.close('all')
        plt.clf()
        print("The figure is saved to {0}".format(sFilename_out))
        return
    else:
        plt.show()
        return
