import os, sys
import numpy as np

from datetime import datetime

from mayavi.mlab import *

sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from eslib.system.define_global_variables import *
from eslib.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def f(x, y):
    sin, cos = np.sin, np.cos
    return sin(x + y) + sin(2 * x - y) + cos(3 * x + 4 * y)
def surface_plot_data_monthly_maya(    aGrid_x, \
    aGrid_y, \
        aData, \
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
                aColor_in = None,\
                aLabel_y_in = None,\
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

    nan_index = np.where(aData == missing_value)
    aData[nan_index] = np.nan
    good_index = np.where(  ~np.isnan(aData))   


    if dMax_z_in is not None:        
        dMax_z = dMax_z_in
    else:         
        dMax_z = np.nanmax(aData) 
    if dMin_z_in is not None:        
        dMin_z = dMin_z_in
    else:
        dMin_z = np.nanmin(aData) #if it has negative value, change here  
    if (dMax_z <= dMin_z ):
        return

    
    
        

    x, y = np.mgrid[-7.:7.05:0.1, -5.:5.05:0.05]
    s = surf(x, y, f)
    #s = surf(aGrid_x, aGrid_y, aData)
    #cs = contour_surf(x, y, f, contour_z=0)
    return s

   

   
