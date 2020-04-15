import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
from scipy.stats import gaussian_kde

sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from eslib.visual.plot.calculate_ticks_space import calculate_ticks_space

def scatter_plot_data_density(aData_x, aData_y,\
    sFilename_out, \
    iSize_x_in = None, \
    iSize_y_in = None,  \
    iDPI_in = None ,\
        dSpace_x_in = None,\
            dSpace_y_in = None,\
    sLabel_x_in =None,\
    sLabel_y_in = None , \
    sLabel_legend_in = None,\
    sTitle_in = None):

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
    if sLabel_x_in is not None:        
        sLabel_X = sLabel_x_in
    else:        
        sLabel_X = ''

    if sLabel_y_in is not None:        
        sLabel_Y = sLabel_y_in
    else:        
        sLabel_Y = ''
    if sLabel_legend_in is not None:        
        sLabel_legend = sLabel_legend_in
    else:        
        sLabel_legend = ''
    if sTitle_in is not None:        
        sTitle = sTitle_in
    else:        
        sTitle = ''

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x )   
    fig.set_figheight( iSize_y )

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005
    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
              
    #ax_scatter = fig.add_axes([0.1, 0.5, 0.8, 0.4] )  
    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    nPoint = len(aData_x)
    y_min = np.nanmin(aData_y) * 0.8#if it has negative value, change here   
    y_max = np.nanmax(aData_y) * 1.2 
        
    x = aData_x
    y = aData_y
    
    df=pd.DataFrame({'x': x, 'y': y })
   
    sns.kdeplot(df.x, df.y, cmap="Reds", shade=True, ax = ax_scatter)
    
    ax_scatter.axis('on')          
    ax_scatter.grid(which='major', color='grey', linestyle='--', axis='y') 
       
    dRatio = 1.0
    ax_scatter.set_aspect(dRatio)  #this one set the y / x ratio
    
    ax_scatter.tick_params(axis="x", labelsize=13) 
    
    ax_scatter.tick_params(axis="y", labelsize=13)
    
    ax_scatter.set_xmargin(0.05)
    ax_scatter.set_ymargin(0.15)
    
    ax_scatter.set_xlabel(sLabel_X,fontsize=12)
    ax_scatter.set_ylabel(sLabel_Y,fontsize=12)
    ax_scatter.set_title( sTitle, loc='center', fontsize=15)
    # round to nearest years...
    
    if y_max < 1000 and y_max > 0.001:
        ax_scatter.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    else: 
        ax_scatter.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
    
    ax_scatter.yaxis.set_major_locator(ticker.MultipleLocator(dSpace))
    y_max = dSpace * 6
    ax_scatter.set_ylim( 0, 75)
    ax_scatter.set_xlim( 0, max(x) )

    dRatio = 75/(max(x)-0.0)
    dRatio = (float(iSize_y)/iSize_x) / ( (75-0 )/ ( max(x)-0.0 ) )
    ax_scatter.set_aspect(dRatio)  #this one set the y / x ratio
    
    #ax_scatter.legend(bbox_to_anchor=(1.0,1.0), loc="upper right",fontsize=12)

    # now determine nice limits by hand:
    
    #binwidthx = 0.01
    #binwidthy = 1
    #limx = np.ceil(np.abs([x]).max() / binwidthx) * binwidthx
    #limy = np.ceil(np.abs([y]).max() / binwidthy) * binwidthy
    #ax_scatter.set_xlim((-lim, lim))
    #ax_scatter.set_ylim((-lim, lim))
    #binsx = np.arange(0, limx + binwidthx, binwidthx)
    #binsy = np.arange(0, limy + binwidthy, binwidthy)
    #ax_histx.hist(x, bins=binsx)
    #ax_histy.hist(y, bins=binsy, orientation='horizontal')
    #ax_histx.set_xlim(ax_scatter.get_xlim())
    #ax_histy.set_ylim(ax_scatter.get_ylim())

    plt.savefig(sFilename_out, bbox_inches='tight')
                       
    plt.close('all')
    print('finished plotting')
