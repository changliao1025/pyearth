
import numpy as np
import matplotlib.pyplot as plt
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.create_line_style import create_line_style


def cdf_plot_multiple_data(aData,
                           sFilename_out,
                           iSize_x_in=None,
                           iSize_y_in=None,
                           iDPI_in=None,
                           dMin_x_in=None,
                           dMax_x_in=None,
                           dSpace_x_in=None,
                           sLabel_x_in=None,
                           sLabel_y_in=None,
                           sTitle_in=None,
                           aLabel_legend_in=None):
    """
    Draw a histogram for single dataset
    """
    nData = len(aData)

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 12
    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 9
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = np.min(aData)

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = np.max(aData)

    if dSpace_x_in is not None:
        dSpace_x = dSpace_x_in
    else:
        # it may be calculated
        pass

    if sLabel_x_in is not None:
        sLabel_x = sLabel_x_in
    else:
        sLabel_x = ''

    if sLabel_y_in is not None:
        sLabel_y = sLabel_y_in
    else:
        sLabel_y = ''

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = np.full(nData, '', dtype=str)

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.15, 0.7
    bottom, height = 0.1, 0.7
    spacing = 0.005
    rect_histogram = [left, bottom, width, height]

    ax_cdf = plt.axes(rect_histogram)
    ax_cdf.tick_params(direction='in', top=True, right=True)

    labels = []
    handles = []
    # setup color
    aColor = create_diverge_rgb_color_hex(nData)
    aLine_style = create_line_style(nData)
    for i in range(0, nData):
        aData_slice = aData[i]
        good_index = np.where((aData_slice >= dMin_x) &
                              (aData_slice <= dMax_x))
        aData_slice = aData_slice[good_index]
        count, bins_count = np.histogram(aData_slice, bins=100)
        pdf = count / sum(count)
        cdf = np.cumsum(pdf)
        ax_cdf.plot(bins_count[1:], cdf, color=aColor[i],
                    label=aLabel_legend[i], linestyle=aLine_style[i])
        labels.append(aLabel_legend[i])
        # handles.append( mpl_patches.Rectangle((0, 0), 1, 1, fc=aColor[i], ec="white", lw=0, alpha=0)  )

    ax_cdf.set_xlim(dMin_x, dMax_x)
    ax_cdf.axis('on')
    ax_cdf.grid(which='major', color='white', linestyle='-', axis='y')

    ax_cdf.set_xlabel(sLabel_x, fontsize=15)
    ax_cdf.set_ylabel(sLabel_y, fontsize=15)

    ax_cdf.legend(labels, loc="lower right", fontsize=15, framealpha=0.7)

    ax_cdf.set_title(sTitle)

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
