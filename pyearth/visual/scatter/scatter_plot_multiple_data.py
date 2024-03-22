
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyearth.visual.scatter.scatter_lowess import scatter_lowess
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex


def scatter_plot_multiple_data(aData_x,
                               aData_y,
                               sFilename_out,
                               iFlag_miniplot_in=None,
                               iFlag_scientific_notation_x_in=None,
                               iFlag_scientific_notation_y_in=None,
                               iSize_x_in=None,
                               iSize_y_in=None,
                               iDPI_in=None,
                               iFlag_log_x_in=None,
                               iFlag_log_y_in=None,
                               dMin_x_in=None,
                               dMax_x_in=None,
                               dMin_y_in=None,
                               dMax_y_in=None,
                               dSpace_x_in=None,
                               dSpace_y_in=None,
                               sFormat_x_in=None,
                               sFormat_y_in=None,
                               sLabel_x_in=None,
                               sLabel_y_in=None,
                               aLabel_point_in=None,
                               aColor_in=None,
                               aMarker_in=None,
                               aLabel_legend_in=None,
                               aSize_in=None,
                               sTitle_in=None,
                               sFont_in=None):
    # number of dataset

    # aData_in = np.array(aData_y)
    try:
        from scipy.stats import gaussian_kde
        from scipy import stats
    except ImportError as e:
        raise ImportError(
            "The package 'scipy' is required for this function to run.") from e

    nData = len(aData_y)

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 12

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 12

    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_miniplot_in is not None:
        iFlag_miniplot = iFlag_miniplot_in
    else:
        iFlag_miniplot = 0

    if iFlag_scientific_notation_x_in is not None:
        iFlag_scientific_notation_x = iFlag_scientific_notation_x_in
    else:
        iFlag_scientific_notation_x = 0

    if iFlag_scientific_notation_y_in is not None:
        iFlag_scientific_notation_y = iFlag_scientific_notation_y_in
    else:
        iFlag_scientific_notation_y = 0

    if iFlag_log_x_in is not None:
        iFlag_log_x = iFlag_log_x_in
    else:
        iFlag_log_x = 0

    if iFlag_log_y_in is not None:
        iFlag_log_y = iFlag_log_y_in
    else:
        iFlag_log_y = 0

    if sLabel_x_in is not None:
        sLabel_X = sLabel_x_in
    else:
        sLabel_X = ''

    if sLabel_y_in is not None:
        sLabel_Y = sLabel_y_in
    else:
        sLabel_Y = ''

    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = ''

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams["font.family"] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    fig = plt.figure(dpi=iDPI)
    fig.set_figwidth(iSize_x)
    fig.set_figheight(iSize_y)

    left, width = 0.1, 0.75
    bottom, height = 0.1, 0.75
    spacing = 0.005
    dY_mini = 0.55
    dX_mini = 0.12
    width_mini = 0.28
    rect_scatter = [left, bottom, width, height]
    rect_scatter_mini = [dY_mini, dX_mini, width_mini, width_mini]
    rect_histx = [left, bottom + height + spacing, width, 0.15]
    rect_histy = [left + width + spacing, bottom, 0.15, height]

    # sns.regplot(x, y, lowess=True)
    # ax_scatter = sns.regplot(x=aData_x, y=aData_y, marker="+", lowess=True)
    ax_scatter_full = plt.axes(rect_scatter)
    if iFlag_miniplot == 1:
        ax_scatter_mini = plt.axes(rect_scatter_mini)
        ax_scatter_all = [ax_scatter_full, ax_scatter_mini]
    else:
        ax_scatter_all = [ax_scatter_full]

    for i in range(nData):
        du = np.array(aData_x[i])
        dv = np.array(aData_y[i])
        du = du.flatten()
        dv = dv.flatten()
        if (i == 0):
            dummyx = du
            dummyy = dv
        else:
            dummyx = np.concatenate([dummyx, du])
            dummyy = np.concatenate([dummyy, dv])

    x_min = np.min([dummyx, dummyy])
    x_max = np.max([dummyx, dummyy])
    y_min = np.min([dummyx, dummyy])
    y_max = np.max([dummyx, dummyy])

    if aColor_in is None:
        # need a better solution here
        if (nData >= 3):
            if (nData <= 10):
                aColor = create_diverge_rgb_color_hex(nData)
            else:
                # we will use both symbol and color
                aColor = create_diverge_rgb_color_hex(5)
                nMarker = np.ceil(nData/5)
                aMarker = ['o', '+', 'x']

                pass
        else:
            if nData == 2:
                aColor = ['red', 'blue']
            else:
                aColor = ['red']
    else:
        aColor = aColor_in

    if aMarker_in is None:
        pass
    else:
        aMarker = aMarker_in

    if aSize_in is None:
        a = mpl.rcParams['lines.markersize'] ** 2
        aSize = np.full(nData, a, dtype=float)
    else:
        aSize = aSize_in

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = x_min

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = x_max

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = y_min

    if dMax_y_in is not None:
        dMax_y = dMax_y_in
    else:
        dMax_y = y_max

    if dSpace_x_in is not None:
        dSpace_x = dSpace_x_in
    else:
        dSpace_x = (dMax_x - dMin_x) / 4.0

    if dSpace_y_in is not None:
        dSpace_y = dSpace_y_in
    else:
        dSpace_y = (dMax_y - dMin_y) / 4.0

    # for ax_scatter in ax_scatter_all:
    dRatio = (float(iSize_y)/iSize_x) / ((dMax_y-dMin_y) / (dMax_x-dMin_x))
    for iax in range(len(ax_scatter_all)):
        ax_scatter = ax_scatter_all[iax]
        ax_scatter.tick_params(direction='in', top=True, right=True)

        if iax == 0:
            aLegend_artist = []
            aLabel = []
        else:
            pass

        for i in range(nData):
            x = aData_x[i]
            y = aData_y[i]

            # treatment for nan data
            x = np.array(x)
            y = np.array(y)
            x = x.flatten()
            y = y.flatten()
            a = np.logical_and(~np.isnan(x), ~np.isnan(y))
            x = x[a]
            y = y[a]

            sc = ax_scatter.scatter(
                x, y, s=aSize[i], alpha=0.5, color=aColor[i], marker=aMarker[i])
            if iax == 0:
                aLegend_artist.append(sc)
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    x, y)
                sR = "slope:" + \
                    "{:.2f}".format(slope) + r"; $r^2$:" + \
                    "{:.2f}".format(r_value**2)
                aLabel.append(aLabel_legend[i] + ': ' + sR)

                if i == 0 and aLabel_point_in is not None:

                    for j in range(len(x)):
                        ax_scatter.annotate(
                            aLabel_point_in[j], (x[j], y[j]), fontsize=16)
                    pass
            # ax_scatter.set_facecolor('silver')

        if iax == 0:

            ax_scatter.axis('on')
            ax_scatter.grid(which='major', color='grey',
                            linestyle='--', axis='y')

            ax_scatter.tick_params(axis="x", labelsize=13)
            ax_scatter.tick_params(axis="y", labelsize=13)

            ax_scatter.set_xmargin(0.05)
            ax_scatter.set_ymargin(0.15)

            ax_scatter.set_xlabel(sLabel_X, fontsize=16)
            ax_scatter.set_ylabel(sLabel_Y, fontsize=16)
            ax_scatter.set_title(sTitle, loc='center', fontsize=20)
            # round to nearest years...

            if sFormat_x_in is not None:
                sFormat_x = sFormat_x_in
                ax_scatter.xaxis.set_major_formatter(
                    mpl.ticker.FormatStrFormatter(sFormat_x))

                # ax_scatter.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1e'))'%.1f'

            if sFormat_y_in is not None:
                sFormat_y = sFormat_y_in
                ax_scatter.yaxis.set_major_formatter(
                    mpl.ticker.FormatStrFormatter(sFormat_y))

                # ax_scatter.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))'%.1f'

            ax_scatter.tick_params(axis='y', pad=8)
            ax_scatter.set_xlim(dMin_x, dMax_x * 1.05)
            ax_scatter.set_ylim(dMin_x, dMax_x * 1.05)
        else:
            ax_scatter.set_xlim(dMin_x, dMax_x * 0.1)
            ax_scatter.set_ylim(dMin_x, dMax_x * 0.1)
            pass

        ax_scatter.xaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='upper', nbins=5))

        if iFlag_log_x == 1:
            aLabel_x = []
            for i in np.arange(dMin_x, dMax_x + 1, 1):
                sTicklabel = r'$10^{{{}}}$'.format(int(i))
                aLabel_x.append(sTicklabel)
                pass

            ax_scatter.set_xticks(np.arange(dMin_x, dMax_x + 1, 1))
            ax_scatter.set_xticklabels(aLabel_x)
        else:
            if iFlag_scientific_notation_x == 1:
                formatter = mpl.ticker.ScalarFormatter(useMathText=True)
                formatter.set_scientific(True)
                # you might need to change here
                formatter.set_powerlimits((-1, 1))
                ax_scatter.xaxis.set_major_formatter(formatter)
            else:
                pass

        # ax_scatter.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base = dSpace_y))
        ax_scatter.yaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='upper', nbins=5))

        if iFlag_log_y == 1:
            aLabel_y = []
            for i in np.arange(dMin_y, dMax_y + 1, 1):
                sTicklabel = r'$10^{{{}}}$'.format(int(i))
                aLabel_y.append(sTicklabel)
                pass

            ax_scatter.set_yticks(np.arange(dMin_y, dMax_y + 1, 1))
            ax_scatter.set_yticklabels(aLabel_y)
            pass
        else:
            if iFlag_scientific_notation_y == 1:
                formatter = mpl.ticker.ScalarFormatter(useMathText=True)
                formatter.set_scientific(True)
                # you might need to change here
                formatter.set_powerlimits((-1, 1))
                ax_scatter.yaxis.set_major_formatter(formatter)

            pass

        ax_scatter.set_aspect(dRatio)  # this one set the y / x ratio

        line = mpl.lines.Line2D([0, 1], [0, 1], color='black', linestyle='dashed')
        transform = ax_scatter.transAxes
        line.set_transform(transform)
        ax_scatter.add_line(line)

        iFlag_lowess = 0
        if (iFlag_lowess == 1):

            y_sm, y_std, order = scatter_lowess(aData_x, aData_y, f=1./3.)
            ax_scatter.plot(x[order], y_sm[order], color='tomato')
            ax_scatter.fill_between(x[order],
                                    y_sm[order] - 1.96*y_std[order],
                                    y_sm[order] + 1.96*y_std[order],
                                    alpha=0.3)

            sLabel_legend_lowess2 = 'LOWESS uncertainty'

            # labels.append(sLabel_legend_lowess2)

        if iax == 0:

            ax_scatter.legend(aLegend_artist, aLabel,
                              loc="upper left", fontsize=16)

        ax_scatter.tick_params(which='both',  # Options for both major and minor ticks
                               top='off',  # turn off top ticks
                               left='off',  # turn off left ticks
                               right='off',  # turn off right ticks
                               bottom='off')  # turn off bottom ticks

        if iax == 0:
            for i in range(nData):
                x = aData_x[i].flatten()
                # treatment for nan data
                a = np.where(~np.isnan(x))
                x = x[a]

                try:
                    if np.max(x) > np.min(x):
                        density = gaussian_kde(x)
                    else:
                        return
                except ValueError:
                    pass

                xx = np.linspace(dMin_x, dMax_x, 1000)
                yy = density(xx)

                # y margin
                y = aData_y[i].flatten()
                a = np.where(~np.isnan(y))
                y = y[a]

                density = gaussian_kde(y)
                xx = np.linspace(dMin_y, dMax_y, 1000)
                yy = density(xx)
                xx, yy = yy, xx

        if iax == 0:

            if iFlag_miniplot == 1:
                # horizontal
                line = mpl.lines.Line2D(
                    [0.0, 0.1], [0.1, 0.1], color='black', linestyle='dotted')
                transform = ax_scatter.transAxes
                line.set_transform(transform)
                ax_scatter.add_line(line)

                # vertical
                line = mpl.lines.Line2D(
                    [0.1, 0.1], [0.1, 0.0], color='black', linestyle='dotted')
                transform = ax_scatter.transAxes
                line.set_transform(transform)
                ax_scatter.add_line(line)

                line = mpl.lines.Line2D([0.1, (dY_mini-left)/width], [
                                     0.1, (dX_mini-bottom + width_mini)/height], color='black', linestyle='dotted')
                transform = ax_scatter.transAxes
                line.set_transform(transform)
                ax_scatter.add_line(line)

                line = mpl.lines.Line2D([0.1,  (dY_mini-left)/width], [0.0,
                                     (dX_mini-bottom)/height], color='black', linestyle='dotted')
                transform = ax_scatter.transAxes
                line.set_transform(transform)
                ax_scatter.add_line(line)

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    print('finished plotting')
