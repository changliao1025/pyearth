import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.choose_n_color import polylinear_gradient, rand_hex_color


def histogram_plot(aData_all,
                   sFilename_out,
                   iSize_x_in=None,
                   iSize_y_in=None,
                   ncolumn_in=None,
                   iFlag_scientific_notation_in=None,
                   iFlag_log_in=None,
                   aColor_in=None,
                   iDPI_in=None,
                   dMin_x_in=None,
                   dMax_x_in=None,
                   dSpace_x_in=None,
                   sLabel_x_in=None,
                   sLabel_y_in=None,
                   sFormat_x_in=None,
                   aLocation_legend_in=None,
                   sLocation_legend_in=None,
                   sTitle_in=None,
                   aLabel_legend_in=None):
    """
    Draw a histogram for single/multiple dataset
    """

    try:
        import scipy
    except ImportError as e:
        raise ImportError(
            "The package 'scipy' is required for this function to run.") from e

    aData_all = np.array(aData_all)
    pShape = aData_all.shape

    nData = pShape[0]

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

    if iFlag_scientific_notation_in is not None:
        iFlag_scientific_notation = iFlag_scientific_notation_in
    else:
        iFlag_scientific_notation = 0

    if iFlag_log_in is not None:
        iFlag_log = iFlag_log_in
    else:
        iFlag_log = 0

    if sFormat_x_in is not None:
        iFlag_format_x = 1
        sFormat_x = sFormat_x_in
    else:
        iFlag_format_x = 0

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = np.min(aData_all)

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = np.max(aData_all)

    if dSpace_x_in is not None:
        iFlag_space_x = 0
        dSpace_x = dSpace_x_in
    else:
        iFlag_space_x = 1
        dSpace_x = (dMax_x - dMin_x) / 20
        pass

    if sLocation_legend_in is not None:
        sLocation_legend = sLocation_legend_in
    else:
        sLocation_legend = "upper right"

    if aLocation_legend_in is not None:
        aLocation_legend = aLocation_legend_in
    else:
        aLocation_legend = (1.0, 1.0)

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
        aLabel_legend = np.full(nData, '')

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
                aColor = ['lightblue']

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.15, 0.7
    bottom, height = 0.1, 0.85
    spacing = 0.005
    rect_histogram = [left, bottom, width, height]

    ax_histo = plt.axes(rect_histogram)
    ax_histo.tick_params(direction='in', top=True, right=True)

    aLegend_artist = []
    aLabel = []

    iFlag_method = 2
    if iFlag_method == 1:  # this method possible won't work for multi-case
        for i in np.arange(1, nData+1):
            aData = aData_all[i-1]
            bad_index = np.where(aData < dMin_x_in)
            aData[bad_index] = dMin_x_in
            bad_index = np.where(aData > dMax_x_in)
            aData[bad_index] = dMax_x_in

            if iFlag_log == 1:
                if dSpace_x >= 1:
                    dSpace_x = int(dSpace_x)
                else:
                    pass

            if i == 0:
                N, bins, hisp = ax_histo.hist(aData, int((dMax_x-dMin_x)/dSpace_x),
                                              color=aColor[i-1],  label=aLabel_legend[i-1])
            else:
                N, bins, hisp = ax_histo.hist(aData, int((dMax_x-dMin_x)/dSpace_x),
                                              color=aColor[i-1],  label=aLabel_legend[i-1], alpha=.5)

            aLegend_artist.append(hisp)
            aLabel.append(aLabel_legend[i-1])
    else:

        bad_index = np.where(aData_all < dMin_x)
        aData_all[bad_index] = dMin_x
        bad_index = np.where(aData_all > dMax_x)
        aData_all[bad_index] = dMax_x

        if iFlag_log == 1:
            if dSpace_x >= 1:
                dSpace_x = int(dSpace_x)
            else:
                pass

        aData_all_t = np.transpose(aData_all)
        N, bins, hisp = ax_histo.hist(aData_all_t, int((dMax_x-dMin_x)/dSpace_x), density=True,
                                      color=aColor, label=aLabel_legend)

        # add density?
        for i in np.arange(1, nData+1):
            aData = aData_all[i-1]
            density = scipy.stats.gaussian_kde(aData)
            good_index = np.where((aData >= dMin_x) & (aData <= dMax_x))
            aData = aData[good_index]
            xx = np.linspace(dMin_x, dMax_x, 1000)
            yy = density(xx)
            ax_histo.plot(xx, yy, color=aColor[i-1], linestyle='--')

    ax_histo.set_xlabel(sLabel_x, fontsize=13)
    ax_histo.set_ylabel(sLabel_y, fontsize=13)
    ax_histo.set_xlim(dMin_x, dMax_x)
    ax_histo.axis('on')

    if iFlag_log == 1:
        aLabel_x = list()
        xtickslocs = ax_histo.get_xticks().tolist()
        for i in np.arange(0, len(xtickslocs), 1):
            ii = xtickslocs[i]
            iii = sFormat_x.format(ii)
            sTicklabel = r'$10^{{{}}}$'.format(iii)
            aLabel_x.append(sTicklabel)
            pass

        ax_histo.set_xticks(xtickslocs)
        ax_histo.set_xticklabels(aLabel_x)
        pass
    else:
        if iFlag_scientific_notation == 1:
            formatter = mpl.ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            ax_histo.yaxis.set_major_formatter(formatter)
            pass
        else:
            if (iFlag_format_x == 1):
                sFormat_x_dummy = sFormat_x.replace("{", "{x")
                ax_histo.xaxis.set_major_formatter(
                    mpl.ticker.StrMethodFormatter(sFormat_x_dummy))
                pass

            pass
        pass

    iFlag_grid = 0  # reserved option
    if iFlag_grid == 1:
        ax_histo.grid(which='major', color='white', linestyle='-', axis='y')

    if iFlag_method == 1:
        ax_histo.legend(aLegend_artist, aLabel, bbox_to_anchor=aLocation_legend,
                        loc=sLocation_legend, fontsize=12)
    else:
        ax_histo.legend(bbox_to_anchor=aLocation_legend,
                        loc=sLocation_legend, fontsize=12)

    ax_histo.set_title(sTitle)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
