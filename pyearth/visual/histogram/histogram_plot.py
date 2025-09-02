import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyearth.visual.color.choose_n_color import polylinear_gradient, rand_hex_color
from pyearth.toolbox.data.list.list_alg import list_max, list_min

def histogram_plot(aData_all_in,
                   sFilename_output_in = None,
                   iSize_x_in=None,
                   iSize_y_in=None,
                   ncolumn_in=None,
                   iFlag_scientific_notation_in=None,
                   iFlag_normalize_in=None, # <-- New parameter for normalization
                   iFlag_log_in=None,
                   iFlag_density_in=None,
                   sHisttype_in=None, # <-- New parameter for histtype
                   aColor_in=None,
                   aPlot_order_in=None,
                   iDPI_in=None,
                   dMin_x_in=None,
                   dMax_x_in=None,
                   dSpace_x_in=None,
                   sLabel_x_in=None,
                   sLabel_y_in=None,
                   sFormat_x_in=None,
                   sFont_in=None,
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

    #make a deep copy of the input

    aData_all = copy.deepcopy(aData_all_in)

    nData = len(aData_all)
    #aData_all = np.array(aData_all)
    #pShape = aData_all.shape
    #nData = pShape[0]

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

    if iFlag_normalize_in is not None: # <-- Initialize iFlag_normalize
        iFlag_normalize = iFlag_normalize_in
    else:
        iFlag_normalize = 0 # Default to False (no normalization)

    if iFlag_density_in is not None:
        iFlag_density = iFlag_density_in
    else:
        iFlag_density = 0

    if iFlag_log_in is not None:
        iFlag_log = iFlag_log_in
    else:
        iFlag_log = 0

    if sFormat_x_in is not None:
        iFlag_format_x = 1
        sFormat_x = sFormat_x_in
    else:
        iFlag_format_x = 0

    if sHisttype_in is not None: # <-- Initialize sHisttype
        sHisttype = sHisttype_in
    else:
        sHisttype = 'bar'

    if dMin_x_in is not None:
        dMin_x = dMin_x_in
    else:
        dMin_x = list_min(aData_all)

    if dMax_x_in is not None:
        dMax_x = dMax_x_in
    else:
        dMax_x = list_max(aData_all)

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

    if sLabel_y_in is not None: # <-- Adjust default y-label based on normalization
        sLabel_y = sLabel_y_in
    else:
        if iFlag_normalize == 1:
            sLabel_y = 'Density'
        else:
            sLabel_y = 'Frequency' # Or 'Counts'

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    if aLabel_legend_in is not None:
        aLabel_legend = aLabel_legend_in
    else:
        aLabel_legend = np.full(nData, '')

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams["font.family"] = sFont
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'

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

    if aPlot_order_in is None:
        aPlot_order = np.arange(nData)
    else:
        # Validate aPlot_order_in
        if not (isinstance(aPlot_order_in, (list, np.ndarray)) and
                len(aPlot_order_in) == nData and
                all(isinstance(idx, (int, np.integer)) for idx in aPlot_order_in) and
                sorted(list(set(aPlot_order_in))) == list(np.arange(nData))):
            raise ValueError(f"aPlot_order_in must be a list or array of unique indices from 0 to {nData-1}, "
                             f"representing a permutation of datasets. Got: {aPlot_order_in}")
        aPlot_order = np.array(aPlot_order_in)

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
    start_alpha = 1.0
    end_alpha = 0.3 if nData > 1 else 1.0
    legend_patch_containers = [None] * nData

    if dSpace_x <= 0:
        if dMax_x == dMin_x :
            num_bins = 1
        else: # dSpace_x is invalid, but range is valid
            num_bins = 20 # Fallback if dSpace_x is bad (original implies 20 if dSpace_x_in is None)
    else:
        calculated_bins = (dMax_x - dMin_x) / dSpace_x
        num_bins = max(1, int(calculated_bins))



    # Main histogram plotting loop, following plot_order
    for plot_sequence_idx, original_idx in enumerate(aPlot_order):
        # Make a copy of the data for this histogram to allow clipping without affecting aData_all
        aData_for_hist = copy.deepcopy(aData_all[original_idx])

        # Clip data for histogram display range
        aData_for_hist[aData_for_hist < dMin_x] = dMin_x
        aData_for_hist[aData_for_hist > dMax_x] = dMax_x

        # Calculate current alpha based on plotting sequence
        current_alpha = start_alpha
        if nData > 1:
            progression = plot_sequence_idx / (nData - 1)
            current_alpha = start_alpha - progression * (start_alpha - end_alpha)
            # Clamp alpha to be within [end_alpha, start_alpha]
            current_alpha = max(end_alpha, min(start_alpha, current_alpha))

        if iFlag_normalize == 1: # <-- Check normalization flag
            N, bins, hisp_patches = ax_histo.hist(aData_for_hist, bins=num_bins,
                                                  color=aColor[original_idx],
                                                  label=aLabel_legend[original_idx],
                                                  alpha=current_alpha,
                                                  density=True,
                                                  histtype=sHisttype) # Normalize to density
        else:
            N, bins, hisp_patches = ax_histo.hist(aData_for_hist, bins=num_bins,
                                                  color=aColor[original_idx],
                                                  label=aLabel_legend[original_idx],
                                                  alpha=current_alpha,
                                                  histtype=sHisttype) # Original call (frequency)
        legend_patch_containers[original_idx] = hisp_patches

    # add density?
    # The density plot loop should also respect aPlot_order for consistent layering
    # and use original unclipped data for KDE fitting.
    # Assuming aPlot_order is correctly used here as well (from previous discussions)
    if iFlag_density == 1:
        for plot_sequence_idx_density, original_idx_density in enumerate(aPlot_order):
            aData_kde = aData_all[original_idx_density] # Use original, unclipped data for KDE
    
            # Add a check for empty or single-value data to prevent KDE errors
            if len(aData_kde) == 0 or (len(aData_kde) > 0 and np.all(aData_kde == aData_kde[0])):
                # print(f"Skipping KDE for dataset {original_idx_density} due to insufficient data or no variance.")
                continue # Skip to the next dataset for KDE
            
            try:
                density = scipy.stats.gaussian_kde(aData_kde)
                xx = np.linspace(dMin_x, dMax_x, 1000)
                yy = density(xx)
                # Determine alpha for density plot (e.g., fixed, or match histogram layer)
                density_plot_alpha = 0.7 # Example: fixed alpha for density lines
                # If you want density alpha to match histogram alpha:
                # density_plot_alpha = start_alpha
                # if nData > 1:
                #     progression_density = plot_sequence_idx_density / (nData - 1)
                #     density_plot_alpha = start_alpha - progression_density * (start_alpha - end_alpha)
                #     density_plot_alpha = max(end_alpha, min(start_alpha, density_plot_alpha))
    
                ax_histo.plot(xx, yy, color=aColor[original_idx_density], linestyle='--', alpha=density_plot_alpha)
            except Exception as e:
                # print(f"Could not compute KDE for dataset {original_idx_density}: {e}")
                pass

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


    #ax_histo.legend(aLegend_artist, aLabel, bbox_to_anchor=aLocation_legend,
    #                    loc=sLocation_legend, fontsize=12)

    final_legend_artists = []
    final_legend_labels = []
    for idx in range(nData): # Iterate in original data order
        if legend_patch_containers[idx] is not None: # Check if this dataset was plotted
            final_legend_artists.append(legend_patch_containers[idx])
            final_legend_labels.append(aLabel_legend[idx])

    if final_legend_artists: # Only show legend if there are items
        ax_histo.legend(final_legend_artists, final_legend_labels, bbox_to_anchor=aLocation_legend,
                            loc=sLocation_legend, fontsize=12)


    ax_histo.set_title(sTitle)

    if sFilename_output_in is None:
        plt.show()
    else:
        #remove it if exists
        if os.path.exists(sFilename_output_in):
            os.remove(sFilename_output_in)

        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == '.png':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            if sExtension == '.pdf':
                plt.savefig(sFilename_out, bbox_inches='tight')
            else:
                plt.savefig(sFilename_out, bbox_inches='tight', format='ps')
        plt.close('all')
        plt.clf()
