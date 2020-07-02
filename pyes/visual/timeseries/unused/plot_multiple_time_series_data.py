import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import AxesGrid

sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from eslib.visual.plot.calculate_ticks_space import calculate_ticks_space
from eslib.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def plot_multiple_time_series_data(aTime, n, \
    aData, sFilename_out, \
    iSize_X_in = None, iSize_Y_in = None,  \
    iDPI_in = None ,\
    sLabel_Y_in = None , aLabel_legend_in = None,\
    sTitle_in = None):

    if iSize_X_in is not None:        
        iSize_X = iSize_X_in
    else:       
        iSize_X = 12
    if iSize_Y_in is not None:        
        iSize_Y = iSize_Y_in
    else:       
        iSize_Y = 9
    if iDPI_in is not None:        
        iDPI = iDPI_in
    else:       
        iDPI = 300
    if sLabel_Y_in is not None:        
        sLabel_Y = sLabel_Y_in
    else:        
        sLabel_Y = ''
    if aLabel_legend_in is not None:        
        aLabel_legend = aLabel_legend_in
    else:        
        aLabel_legend = np.full(n,'')
    if sTitle_in is not None:        
        sTitle = sTitle_in
    else:        
        sTitle = ''

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_X)   
    fig.set_figheight( iSize_Y)
              
    pAxGrid = AxesGrid(fig, 111,
            nrows_ncols=(1,1),
            axes_pad=0.6,    
            label_mode='')  # note the empty label_mode

    pYear = mdates.YearLocator()   # every year
    pMonth = mdates.MonthLocator()  # every month
    sYear_format = mdates.DateFormatter('%Y')

    nstress = len(aTime)
    y_min = np.nanmin(aData) #if it has negative value, change here   
    y_max = np.nanmax(aData) * 1.2 

    for i, ax in enumerate(pAxGrid):      

        aColor = create_diverge_rgb_color_hex(n)
        for i in range(n):
            x1 = aTime
            y1 = aData[i]
            ax.plot( x1, y1, \
                    color = aColor[i], linestyle = '--' , \
                    label= aLabel_legend[i])
        
        #now let's set the axis
        ax.axis('on')          

        ax.xaxis.set_major_locator(pYear)        
        ax.xaxis.set_minor_locator(pMonth)
        ax.xaxis.set_major_formatter(sYear_format)

        ax.tick_params(axis="x", labelsize=13)
        ax.tick_params(axis="y", labelsize=13)
        ax.set_xmargin(0.05)
        ax.set_ymargin(0.15)
        
        ax.set_xlabel('Year',fontsize=12)
        ax.set_ylabel(sLabel_Y,fontsize=12)

        ax.set_title( sTitle, loc='center', fontsize=15)

        # round to nearest years...
        x_min = np.datetime64(aTime[0], 'Y')
        x_max = np.datetime64(aTime[nstress-1], 'Y') + np.timedelta64(1, 'Y')

        ax.set_xlim(x_min, x_max)
        if y_max < 1000 and y_max > 0.001:
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        else: 
            ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
        dSpace = calculate_ticks_space(aData)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(dSpace))
        y_max = dSpace * 6
        ax.set_ylim( y_min, y_max)
        ax.grid(which='major', color='grey', linestyle='--', axis='y') 
        #ax.grid(which='minor', color='#CCCCCC', linestyle=':') #only y     axis grid is 

        dRatio = (float(iSize_Y)/iSize_X) / (   (y_max-y_min )/ nstress )
        ax.set_aspect(dRatio)  #this one set the y / x ratio
        
        ax.legend(bbox_to_anchor=(1.0,1.0), loc="upper right",fontsize=12)

    plt.savefig(sFilename_out, bbox_inches='tight')
                       
    plt.close('all')
    print('finished plotting')
