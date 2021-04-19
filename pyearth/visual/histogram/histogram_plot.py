import os, sys
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpl_patches

def histogram_plot(aData, \
    sFilename_out, \
    iSize_x_in = None, \
    iSize_y_in = None, \
    iDPI_in = None, \
    dMin_in = None, \
    dMax_in = None, \
    dMin_x_in = None, \
    dMax_x_in = None, \
    dSpace_x_in = None, \
    dSpace_y_in = None, \
    sLabel_x_in = None, \
    sLabel_y_in = None, \
    sTitle_in = None):

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
    if dMin_in is not None:        
        dMin = dMin_in
    else:       
        dMin = np.min(aData)

    if dMax_in is not None:        
        dMax = dMax_in
    else:       
        dMax = np.max(aData)

    if dMin_x_in is not None:        
        dMin_x = dMin_x_in
    else:       
        dMin_x = dMin

    if dMax_x_in is not None:        
        dMax_x = dMax_x_in
    else:       
        dMax_x = dMax

    if dSpace_x_in is not None:        
        dSpace_x = dSpace_x_in
    else:       
        #it may be calculated 
        pass

    if dSpace_y_in is not None:        
        dSpace_y = dSpace_y_in
    else:       
        #it may be calculated 
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
    
    good_index = np.where( (aData >= dMin) & (aData<= dMax)  )
    
    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )   
    fig.set_figheight( iSize_y )

    left, width = 0.15, 0.7
    bottom, height = 0.1, 0.85
    spacing = 0.005
    rect_histogram = [left, bottom, width, height]
   
    ax_histo = plt.axes(rect_histogram)
    ax_histo.tick_params(direction='in', top=True, right=True)

    aData = aData[good_index]

    ax_histo.hist(aData, int((dMax_x-dMin_x)/dSpace_x))  

    ax_histo.set_xlabel(sLabel_x,fontsize=13 )
    ax_histo.set_ylabel(sLabel_y,fontsize=13 )   
   
  
    ax_histo.set_xlim( dMin_x, dMax_x )
   
    ax_histo.axis('on')   
    ax_histo.grid(which='major', color='white', linestyle='-', axis='y')
    
   
    leg = ax_histo.legend(bbox_to_anchor=(1.0,1.0),loc='upper right', frameon=True)
    
    frame = leg.get_frame()
    frame.set_edgecolor('black')

    ax_histo.set_title(sTitle)
    
    plt.savefig(sFilename_out, bbox_inches='tight')
                       
    plt.close('all')
    print('finished plotting')