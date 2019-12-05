import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import AxesGrid

def y_formatter(x):
            return '$%1.2f' % x

def plot_time_series_data(aTime, aData, sFilename_out, iSize_X_in = None, iSize_Y_in = None,  \
iDPI_in = None , sLabel_Y_in = None , sTitle_in = None):


    if iSize_X_in is not None:
        
        iSize_X = iSize_X_in
    else:
        #by default, this model will run in steady state
        iSize_X = 12
    if iSize_Y_in is not None:
        
        iSize_Y = iSize_Y_in
    else:
        #by default, this model will run in steady state
        iSize_Y = 9

    if iDPI_in is not None:        
        iDPI = iDPI_in
    else:
        #by default, this model will run in steady state
        iDPI = 300
    if sLabel_Y_in is not None:        
        sLabel_Y = sLabel_Y_in
    else:
        #by default, this model will run in steady state
        sLabel_Y = ''
    if sTitle_in is not None:        
        sTitle = sTitle_in
    else:
        #by default, this model will run in steady state
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
    y_min = np.nanmin(aData)
   
    y_max = np.nanmax(aData) * 1.2 


    for i, ax in enumerate(pAxGrid):                   
        ax.axis('on')   
        

        ax.xaxis.set_major_locator(pYear)        
        ax.xaxis.set_minor_locator(pMonth)
        ax.xaxis.set_major_formatter(sYear_format)

        
        ax.set_xmargin(0.05)
        ax.set_ymargin(0.15)
        
        ax.set_ylim( y_min, y_max)
        
        ax.set_xlabel('Year')
        ax.set_ylabel(sLabel_Y)

        ax.set_title( sTitle, loc='center')

        # round to nearest years...
        x_min = np.datetime64(aTime[0], 'Y')
        x_max = np.datetime64(aTime[nstress-1], 'Y') + np.timedelta64(1, 'Y')

        ax.set_xlim(x_min, x_max)
       
        #ax.format_ydata = y_formatter
        x1 = aTime
        y1 = aData
        
        ax.grid(which='major', color='#CCCCCC', linestyle='--', axis='y') 
        #ax.grid(which='minor', color='#CCCCCC', linestyle=':')
        
        dRatio = (float(iSize_Y)/iSize_X) / (   (y_max-y_min )/ nstress )
        ax.set_aspect(dRatio)  #this one set the y / x ratio
        ax.plot( x1, y1, color = 'red', linestyle = '--' , marker="+", markeredgecolor='blue' , label='Observed precipitation')
        ax.legend(bbox_to_anchor=(0.95,0.95), loc="upper right")

    plt.savefig(sFilename_out, bbox_inches='tight')
                       
    plt.close('all')