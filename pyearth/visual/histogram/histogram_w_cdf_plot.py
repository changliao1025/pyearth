
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def histogram_w_cdf_plot(aData,
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
                         sLabel_legend_in=None):
    """
    Draw a histogram for single dataset
    """

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

    if sLabel_legend_in is not None:
        sLabel_legend = sLabel_legend_in
    else:
        sLabel_legend = ''

    good_index = np.where((aData >= dMin_x) & (aData <= dMax_x))

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.15, 0.7
    bottom, height = 0.1, 0.85
    spacing = 0.005
    rect_histogram = [left, bottom, width, height]

    ax_histo = plt.axes(rect_histogram)
    ax_histo.tick_params(direction='in', top=True, right=True)
    aData = aData[good_index]
    ax_histo.hist(aData, int((dMax_x-dMin_x)/dSpace_x), density=True)
    ax_histo.set_xlabel(sLabel_x, fontsize=13)
    ax_histo.set_ylabel(sLabel_y, fontsize=13)
    ax_histo.set_xlim(dMin_x, dMax_x)

    ax_histo.axis('on')
    ax_histo.grid(which='major', color='white', linestyle='-', axis='y')

    handles = [mpl.patches.Rectangle(
        (0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 1

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append(sLabel_legend)
    # create the legend, supressing the blank space of the empty line symbol    and the
    # padding between symbol and label by setting handlelenght and  handletextpad

    ax_histo.legend(handles, labels, loc="upper right", fontsize=12,
                    fancybox=True, framealpha=0.7,
                    handlelength=0, handletextpad=0)

    ax_cdf = ax_histo.twinx()
    count, bins_count = np.histogram(aData, bins=100)
    # finding the PDF of the histogram using count values
    pdf = count / sum(count)
    cdf = np.cumsum(pdf)
    ax_cdf.plot(bins_count[1:], cdf)

    ax_histo.set_title(sTitle)
    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
