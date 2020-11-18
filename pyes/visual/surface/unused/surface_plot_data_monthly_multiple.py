import os, sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from datetime import datetime
from mpl_toolkits import mplot3d
from matplotlib.collections import PolyCollection
from matplotlib import cm


sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from pyes.system.define_global_variables import *
from pyes.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def surface_plot_data_monthly_multiple(aData_x, \
    aData_y, \
    aData_all, \
    sFilename_out,\
    iDPI_in = None,\
    iFlag_trend_in = None, \
    iReverse_z_in = None, \
    iSize_x_in = None, \
    iSize_y_in = None, \
    dMin_x_in=None,\
    dMax_x_in=None,\
    dMax_z_in =None, \
    dMin_z_in = None, \
    dSpace_x_in =None,\
    dSpace_z_in=None,\
    sMarker_in =None,\
    sLabel_x_in = None, \
    sLabel_y_in = None, \
    sLabel_z_in = None, \
    sLabel_legend_in = None,\
    sTitle_in = None):

    if iDPI_in is not None:        
        iDPI = iDPI_in
    else:       
        iDPI = 300
    
    if iFlag_trend_in is not None:
        iFlag_trend = 1
    else:
        iFlag_trend = 0
        
    if iReverse_z_in is not None:
        iReverse_z = 1
    else:
        iReverse_z = 0

    if iSize_x_in is not None:        
        iSize_x = iSize_x_in
    else:       
        iSize_x = 12
    if iSize_y_in is not None:        
        iSize_y = iSize_y_in
    else:       
        iSize_y = 9
    if dSpace_x_in is not None:        
        dSpace_x = dSpace_x_in
    else:       
        dSpace_x = 1
    if dSpace_z_in is not None:        
        dSpace_z = dSpace_z_in
    else:       
        dSpace_z = 1

    
    if sLabel_x_in is not None:        
        sLabel_x = sLabel_x_in
    else:        
        sLabel_x = ''

    if sLabel_y_in is not None:        
        sLabel_y = sLabel_y_in
    else:        
        sLabel_y = ''
    if sLabel_z_in is not None:        
        sLabel_z = sLabel_z_in
    else:        
        sLabel_z = ''
    if sLabel_legend_in is not None:        
        sLabel_legend = sLabel_legend_in
    else:        
        sLabel_legend = ''
    if sTitle_in is not None:        
        sTitle = sTitle_in
    else:        
        sTitle = ''

    if sMarker_in is not None:        
        sMarker = sMarker_in
    else:        
        sMarker = '+'

  

    nslice = len(aData_all)
    

    if dMax_z_in is not None:        
        dMax_z = dMax_z_in
    else:         
        dMax_z = np.nanmax(aData_all) 
    if dMin_z_in is not None:        
        dMin_z = dMin_z_in
    else:
        dMin_z = np.nanmin(aData_all) #if it has negative value, change here  
    if (dMax_z <= dMin_z ):
        return

    fig = plt.figure( dpi=iDPI )
    fig.set_figwidth( iSize_x)   
    fig.set_figheight( iSize_y)         

    ax = fig.add_subplot(111, projection='3d')
    #ax.pbaspect = np.array([2.0, 10.0, 0.5])    
    ax.view_init(40,300)
    
    verts=[]
    ys = range(nslice)
    for iSlice in range(2, nslice+1):
        aData = aData_all[iSlice-1]
        nan_index = np.where(aData==missing_value)
        #aData[nan_index] = -9999
        good_index = np.where(  ~np.isnan(aData))   
        
        # Make data.
        #X = np.arange(-5, 5, 0.25)
        #Y = np.arange(-5, 5, 0.25)
        #X, Y = np.meshgrid(X, Y)
        X = aData_x
        Y = aData_y
        
        Z = aData
        surf = ax.plot_surface( X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
        
    ax.axis('on')              

    ax.xaxis._axinfo["grid"]['linewidth'] = 0.
    ax.yaxis._axinfo["grid"]['linewidth'] = 0.
    ax.zaxis._axinfo["grid"]['color'] = "grey"
    ax.zaxis._axinfo["grid"]['linestyle'] = "--"

    

    ax.set_xlim3d( np.min(aData_x), np.max(aData_x) )
    ax.set_xlabel(sLabel_x)
    

    ax.set_ylim3d( np.min(aData_y), np.max(aData_y) )    
    ax.set_ylabel(sLabel_y)

    if (iReverse_z==1):
        ax.set_zlim3d(80, 0)
    else:
        ax.set_zlim3d(np.nanmin(aData), np.nanmax(aData))

    ax.set_zlabel(sLabel_z)

    plt.savefig(sFilename_out, bbox_inches='tight')
                       
    plt.close('all')
    plt.clf()
    #print('finished plotting')
