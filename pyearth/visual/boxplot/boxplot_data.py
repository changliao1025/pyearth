import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

mpl.rcParams['hatch.linewidth'] = 0.5


def boxplot_data(aData_in,
                 aLabel_x_in,
                 aLabel_y_in,
                 sFilename_out,
                 iDPI_in=None,
                 ncolumn_in=None,
                 iSize_x_in=None,
                 iSize_y_in=None,
                 dMax_y_in=None,
                 dMin_y_in=None,
                 dSpace_y_in=None,
                 aMarker_in=None,
                 aColor_in=None,
                 aHatch_in=None,
                 sLabel_y_in=None,
                 aLabel_legend_in=None,
                 aLocation_legend_in=None,
                 sFormat_y_in=None,
                 sLocation_legend_in=None,
                 sTitle_in=None):
    """Draw a boxplot

    Args:
        aData_in (_type_): _description_
        aLabel_x_in (_type_): _description_
        aLabel_y_in (_type_): _description_
        sFilename_out (_type_): _description_
        iDPI_in (_type_, optional): _description_. Defaults to None.
        ncolumn_in (_type_, optional): _description_. Defaults to None.
        iSize_x_in (_type_, optional): _description_. Defaults to None.
        iSize_y_in (_type_, optional): _description_. Defaults to None.
        dMax_y_in (_type_, optional): _description_. Defaults to None.
        dMin_y_in (_type_, optional): _description_. Defaults to None.
        dSpace_y_in (_type_, optional): _description_. Defaults to None.
        aMarker_in (_type_, optional): _description_. Defaults to None.
        aColor_in (_type_, optional): _description_. Defaults to None.
        aHatch_in (_type_, optional): _description_. Defaults to None.
        sLabel_y_in (_type_, optional): _description_. Defaults to None.
        aLabel_legend_in (_type_, optional): _description_. Defaults to None.
        aLocation_legend_in (_type_, optional): _description_. Defaults to None.
        sFormat_y_in (_type_, optional): _description_. Defaults to None.
        sLocation_legend_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
    """
    # if a data does not have a consistent dimension, can we use boxplot?
    # Yes, in this case, a list of list is a better data structure for the data

    nData = len(aData_in)

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
        dMax_y = -9999
        for i in range(nData):
            aData_dummy = aData_in[i]
            for j in range(len(aData_dummy)):
                dMax_y_dummy = np.max(aData_dummy[j])
                dMax_y = np.max([dMax_y, dMax_y_dummy])

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = 9999
        for i in range(nData):
            aData_dummy = aData_in[i]
            for j in range(len(aData_dummy)):
                dMin_y_dummy = np.min(aData_dummy[j])
                dMin_y = np.min([dMin_y, dMin_y_dummy])

    if (dMax_y <= dMin_y):
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
        aLocation_legend = (1.0, 1.0)

    if ncolumn_in is not None:
        ncolumn = ncolumn_in
    else:
        ncolumn = 1

    if aColor_in is not None:
        aColor = aColor_in
    else:
        if (nData >= 3):
            aColor = create_diverge_rgb_color_hex(nData)
        else:
            if nData == 2:
                aColor = ['red', 'blue']
            else:
                aColor = ['red']

    if aHatch_in is not None:
        aHatch = aHatch_in
    else:
        aHatch = np.full(nData, '+')

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])

    dWidth_total = 0.7

    aLegend_artist = []
    aLabel = []

    for i in np.arange(1, nData + 1, 1):

        data1 = aData_in[i-1]

        # the data1 is still considered a list
        nData1 = len(data1)
        # plot group

        # set x tick
        x_center = i-1
        # set width
        dWidth = 0.8 * dWidth_total / nData1

        for j in range(1, nData1+1, 1):
            # plot individual data
            data_raw = data1[j-1]

            q1_raw = np.percentile(data_raw, 25)
            q3_raw = np.percentile(data_raw, 75)
            qm_raw = np.median(data_raw)
            #qm = qm_raw  # np.log10(qm_raw)

            iqr_raw = q3_raw - q1_raw
            #iqr_log = np.log10(iqr_raw)

            q0_raw = np.percentile(data_raw, 0.35)
            q4_raw = np.percentile(data_raw, 99.65)  # q3_raw + 1.5 * iqr_raw

            # set
            x_left = x_center - 0.5 * dWidth_total + (j-1) * dWidth
            x_right = x_left + dWidth
            x_mid = 0.5*(x_left + x_right)

            y_bot = q1_raw
            y_top = q3_raw
            height = y_top-y_bot
            # draw box
            rect = mpl.patches.Rectangle((x_left, y_bot), dWidth, height, linewidth=0.5,
                                     facecolor=aColor[j-1], edgecolor='k')  # , hatch = aHatch[j-1 ]

            if i == 1:
                aLegend_artist.append(rect)
                aLabel.append(aLabel_legend_in[j-1])

            ax.add_patch(rect)
            # draw lower q1
            x0 = [x_mid, x_mid]
            y0 = [q0_raw, q1_raw]
            line, = ax.plot(x0, y0,
                            color='black', linestyle='--', linewidth=1)
            # draw short line
            x0 = [x_mid - 0.2 * dWidth, x_mid + 0.2 * dWidth]
            y0 = [q0_raw, q0_raw]
            line, = ax.plot(x0, y0,
                            color='black', linestyle='-', linewidth=1)

            # draw upper q3
            x0 = [x_mid, x_mid]
            y0 = [q3_raw, q4_raw]
            line, = ax.plot(x0, y0,
                            color='black', linestyle='--', linewidth=1)

            # draw short line
            x0 = [x_mid - 0.2 * dWidth,  x_mid + 0.2 * dWidth]
            y0 = [q4_raw, q4_raw]
            line, = ax.plot(x0, y0,
                            color='black', linestyle='-', linewidth=1)

            # draw median
            x0 = [x_mid - 0.5 * dWidth,  x_mid + 0.5 * dWidth]
            y0 = [qm_raw, qm_raw]
            line, = ax.plot(x0, y0,
                            color='black', linestyle='-', linewidth=1)

            # whiskers outlier
            dummy_index = np.where(data_raw < q0_raw)[0]
            for k in range(len(dummy_index)):
                x0 = [x_mid]
                y0 = [data_raw[dummy_index[k]]]
                ax.plot(x0, y0, marker="o", markersize=3,
                        markeredgecolor="black", markerfacecolor="white")
                pass
            dummy_index = np.where(data_raw > q4_raw)[0]
            for k in range(len(dummy_index)):
                x0 = [x_mid]
                y0 = [data_raw[dummy_index[k]]]
                ax.plot(x0, y0, marker="o", markersize=3,
                        markeredgecolor="black", markerfacecolor="white")
                pass

        pass

    x = np.arange(len(aLabel_x_in))
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(sLabel_y)
    ax.set_title(sTitle)
    ax.set_xticks(x)
    ax.set_xticklabels(aLabel_x_in)

    if (iFlag_format_y == 1):
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(sFormat_y))

    ax.set_xlim(-1, nData)
    ax.set_ylim(dMin_y, dMax_y)
    ax.grid(which='major', color='grey', linestyle='--', axis='y')
    #ax.set_yscale('log')

    ax.legend(aLegend_artist, aLabel, bbox_to_anchor=aLocation_legend,
              loc=sLocation_legend, fontsize=12, ncol=ncolumn)

    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close('all')
    plt.clf()
