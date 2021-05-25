import matplotlib.pyplot as plt
import numpy as np

def plot_time_series_vertical_data(aTime, \
    aVertical,\
        aData,\
            sFilename_out,\
                iDPI_in = None,\
                    iFlag_log_in = None,\
      iReverse_y_in = None    ,\
      iSize_x_in = None, \
                          iSize_y_in = None ,\
                              dMax_x_in = None, \
                          dMin_x_in = None, \
                          dMax_y_in = None, \
                          dMin_y_in = None,\
                              sLabel_y_in=None,\
                                  sLabel_legend_in = None,\
                                             sLocation_legend_in=None,\
                                                 aLocation_legend_in =None):
    # Generate data for the plot
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_log_in is not None:
        iFlag_log = 1
    else:
        iFlag_log = 0

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
        dMax_y = np.nanmax(aVertical) #* 1.2

    if dMin_y_in is not None:
        dMin_y = dMin_y_in
    else:
        dMin_y = np.nanmin(aVertical) #if it has negative value, change here

    if (dMax_y <= dMin_y ):
        return

    if sLabel_y_in is not None:
        sLabel_y = sLabel_y_in
    else:
        sLabel_y = ''

    if sLabel_legend_in is not None:
        iFlag_legend=1
        sLabel_legend = sLabel_legend_in
    else:
        iFlag_legend=0
        sLabel_legend = ''
    if sLocation_legend_in is not None:
        sLocation_legend = sLocation_legend_in
    else:
        sLocation_legend = "upper right"

    if aLocation_legend_in is not None:
        aLocation_legend = aLocation_legend_in
    else:
        aLocation_legend=(1.0,1.0)

    # Generate the plot
    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.5, 0.8, 0.4] )
    if (iReverse_y ==1): #be careful here
        ax.set_ylim( dMax_y, dMin_y )
    else:
        ax.set_ylim( dMin_y, dMax_y )
    cmap = ax.pcolormesh(aTime, aVertical, aData)

    fig.colorbar(cmap)

    #next Y axis
    ax.set_ylabel(sLabel_y,fontsize=12)
    if iFlag_legend ==1 :
        ax.legend(bbox_to_anchor=aLocation_legend, \
              loc=sLocation_legend, \
              fontsize=12)

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    return