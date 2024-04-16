import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex


def barplot_data_with_reference(aData_in,
                                aLabel_x_in,
                                aLabel_y_in,
                                sFilename_out,
                                aLabel_z_in=None,
                                aData_reference_in=None,
                                aLabel_legend_reference_in=None,
                                iDPI_in=None,
                                iFlag_scientific_notation_in=None,
                                ncolumn_in=None,
                                iSize_x_in=None,
                                iSize_y_in=None,
                                dMax_y_in=None,
                                dMin_y_in=None,
                                dSpace_y_in=None,
                                aMarker_in=None,
                                aColor_in=None,
                                aHatch_in=None,
                                sLabel_info_in=None,
                                sLabel_y_in=None,
                                aLinestyle_in=None,
                                aLocation_legend_in=None,
                                sFormat_y_in=None,
                                sLocation_legend_in=None,
                                sTitle_in=None):

    aData_in = np.array(aData_in)
    aData_reference_in = np.array(aData_reference_in)
    pShape = aData_in.shape
    ndim = aData_in.ndim
    if ndim == 2:
        iFlag_sub = 0
        nCat = pShape[0]
        nData = pShape[1]

    else:
        iFlag_sub = 1
        nCat = pShape[0]
        nsub = pShape[1]
        nData = pShape[2]

    if aData_reference_in is not None:
        iFlag_ref = 1
        nData_reference = aData_reference_in.shape[0]

    else:
        iFlag_ref = 0

    if iFlag_scientific_notation_in is not None:
        iFlag_scientific_notation = 1
    else:
        iFlag_scientific_notation = 0

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
        dMin_y = np.nanmin(aData_in)  # if it has negative value, change here

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
            aColor = create_diverge_rgb_color_hex(nData+nData_reference)
        else:
            if nData == 2:
                aColor = ['red', 'blue']
            else:
                aColor = ['red']

    if aHatch_in is not None:
        aHatch = aHatch_in
    else:
        aHatch = np.fill(nData, '+')

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4])
    x = np.arange(len(aLabel_x_in))
    dMin_x = -0.5
    dMax_x = len(aLabel_x_in)-0.5

    ax.set_ylabel(sLabel_y, fontsize=14)
    ax.set_title(sTitle, fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(aLabel_x_in)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)

    total_width = 0.6

    aLegend_artist = []
    aLabel = []
    if ndim == 2:
        width = total_width / (nData)
        for i in range(0, nCat, 1):
            for j in np.arange(0, nData, 1):
                data1 = aData_in[j, i]
                x1 = x[i] - total_width * 0.5 + (j+0.5) * width
                rects = ax.bar(x1, data1, width, label=aLabel_y_in[k], linestyle=aLinestyle_in[k],
                               color=aColor[k], hatch=aHatch[k], edgecolor="k")

                pass
        if iFlag_ref == 1:
            for i in aData_reference_in:
                x0 = [-1, nData-len(aData_reference_in)]
                y0 = [aData_in[i][0], aData_in[i][0]]
                ax.plot(x0, y0,
                        color=aColor[i], linestyle='dashed',
                        marker=aMarker_in[i],
                        label=aLabel_y_in[i])
    else:
        width = total_width / (nData * nsub)
        if iFlag_ref == 1:
            for i in range(nData_reference):
                x0 = [-1, nCat]
                y0 = [aData_reference_in[i], aData_reference_in[i]]
                line, = ax.plot(x0, y0,
                                color=aColor[i+nData], linestyle='dashed',
                                marker=aMarker_in[i],
                                label=aLabel_legend_reference_in[i])
                aLegend_artist.append(line)
                aLabel.append(aLabel_legend_reference_in[i])

        for i in range(0, nCat, 1):
            for j in np.arange(0, nsub, 1):
                for k in np.arange(0, nData, 1):
                    data1 = aData_in[i, j, k]
                    x1 = x[i] - total_width * 0.5 + j*width + k * (width * 2)
                    # print(x1)
                    if j == 0 and i == 0:
                        rects, = ax.bar(x1, data1, width, label=aLabel_y_in[k], linestyle=aLinestyle_in[k],
                                        color=aColor[k], hatch=aHatch[j], edgecolor="k")
                        aLegend_artist.append(rects)
                        aLabel.append(aLabel_y_in[k])

                    else:
                        rects, = ax.bar(x1, data1, width,  linestyle=aLinestyle_in[k],
                                        color=aColor[k], hatch=aHatch[j], edgecolor="k")

        pass

    if iFlag_scientific_notation == 1:
        formatter = mpl.ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        # formatter.set_powerlimits((-1,1)) # you might need to change here
        ax.yaxis.set_major_formatter(formatter)
        # most time, when you use scientific notation, you may not need set the space,
        # but you may still set it using the method below
        pass

    if sLabel_info_in is not None:
        ax.text(0.1, 0.9, sLabel_info_in,
                verticalalignment='center', horizontalalignment='left',
                transform=ax.transAxes,
                color='black', fontsize=13)

    if (iFlag_format_y == 1):
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(sFormat_y))

    ax.set_xlim(dMin_x, dMax_x)
    ax.set_ylim(dMin_y, dMax_y)
    ax.grid(linewidth=1, color='gray', alpha=0.3, linestyle='--')

    if iFlag_sub == 1:
        handles = list()
        # labels = []
        for i in range(nsub):
            # handles.append(
            p = mpl.patches.Rectangle(
                (0, 0), 2, 2, hatch=aHatch[i], facecolor='w', label=aLabel_z_in[i])  # )
            aLegend_artist.append(p)
            aLabel.append(aLabel_z_in[i])

    ax.legend(aLegend_artist, aLabel,      bbox_to_anchor=aLocation_legend,
              loc=sLocation_legend,
              fontsize=14,
              ncol=ncolumn)

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()

    return
