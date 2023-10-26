import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.stats import gaussian_kde

sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
def ridgeplot_data_density(aDict, 
    aData, 
        sFilename_out,
            iSize_x_in = None, 
                              iSize_y_in = None, 
                              iDPI_in = None,
                              iFlag_scientific_notation_x_in=None,
                              iFlag_scientific_notation_y_in=None,
                              iFlag_log_x_in=None,
                              iFlag_log_y_in=None,
                              dMin_x_in = None, 
                              dMax_x_in = None, 
                              dSpace_x_in = None, 
                              sFormat_x_in =None,
                              sLabel_x_in = None, 
                              sLabel_y_in = None, 
                              sLabel_legend_in = None, 
                              sTitle_in = None):
    nData = len(aDict)
    
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

    if iFlag_log_x_in is not None:
        iFlag_log_x = iFlag_log_x_in
    else:
        iFlag_log_x = 0

    if iFlag_log_y_in is not None:
        iFlag_log_y = iFlag_log_y_in
    else:
        iFlag_log_y = 0

    if iFlag_scientific_notation_x_in is not None:
        iFlag_scientific_notation_x = iFlag_scientific_notation_x_in
    else:
        iFlag_scientific_notation_x = 0

    if iFlag_scientific_notation_y_in is not None:
        iFlag_scientific_notation_y = iFlag_scientific_notation_y_in
    else:
        iFlag_scientific_notation_y = 0

    if sLabel_x_in is not None:
        sLabel_x = sLabel_x_in
    else:
        sLabel_x = ''

    if sLabel_y_in is not None:
        sLabel_y = sLabel_y_in
    else:
        sLabel_y = ''

    if sLabel_legend_in is not None:
        sLabel_legend = sLabel_legend_in
    else:
        sLabel_legend = ''

    if sTitle_in is not None:
        sTitle = sTitle_in
    else:
        sTitle = ''

    plt.rcParams["font.family"] = "Times New Roman"
    fig = plt.figure( dpi=iDPI )
    

    # we generate a color palette with Seaborn.color_palette()
    pal = sns.color_palette(palette='coolwarm', n_colors=nData)

    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y ) 
    axgr = AxesGrid(fig, 111, 
                    nrows_ncols=(nData, 1),
                    axes_pad=-0.1,      
                    label_mode='')  # note the empty label_mode

                   

    if dMin_x_in is None:
        for i in range(nData):
            aData_dummy = np.array(aData[i]   )
            dMin_x_dummy = np.min(aData_dummy)          
            if i ==0:
                dMin_x = dMin_x_dummy              
            else:
                dMin_x = np.min([dMin_x, dMin_x_dummy] )
              
    else:
        dMin_x = dMin_x_in
    
    if dMax_x_in is None:
        for i in range(nData):
            aData_dummy = np.array(aData[i]   )           
            dMax_x_dummy = np.max(aData_dummy)
            if i ==0:              
                dMax_x = dMax_x_dummy
            else:               
                dMax_x = np.max([dMax_x, dMax_x_dummy])
    else:
        dMax_x = dMax_x_in

    

    for i, ax in enumerate(axgr):
        aData_dummy = aData[i]        

        density = gaussian_kde(aData_dummy)
        xx = np.linspace(dMin_x, dMax_x,1000)
        yy = density(xx)
        ax.plot(xx,yy, color='w',linewidth=0.5)#
        ax.fill_between(xx, yy, 0,  linewidth=1.5, color = pal[i])
        dMin_y = np.min(yy)
        dMax_y = np.max(yy)
        
        dRatio = (float((iSize_y/nData))/iSize_x) / ( (dMax_y-dMin_y )/ ( dMax_x-dMin_x ) )
        ax.set_aspect(dRatio, 'box')

        ax.spines.left.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.spines.bottom.set_visible(False)        

        ax.get_yaxis().set_visible(False)
        sText = aDict[i+1]
        ax.text(0.85, 0.40, sText, 
        verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, 
            color= pal[i], fontsize=10, fontweight='bold')

        if i < (nData-1):
            
            ax.axes.xaxis.set_visible(False)
        else:  
            sText = sLabel_x
            ax.set_xlabel(sText,  fontsize=15)
            pass
    
   
    

    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    return