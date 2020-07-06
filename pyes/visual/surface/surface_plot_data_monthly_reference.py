import os, sys
import numpy as np

from datetime import datetime
import matplotlib.pyplot as plt
import pyvista as pv


sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from pyes.system.define_global_variables import *
from pyes.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex


def surface_plot_data_monthly_reference(    aGrid_x, \
    aGrid_y, \
    aData, \
    aData_reference,\
    sFilename_out,\
    iDPI_in = None,\
    iFlag_trend_in = None, \
    iReverse_z_in = None, \
    iSize_x_in = None, \
    iSize_y_in = None, \
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
        dMax_z = np.nanmax(aData_reference) 
    if dMin_z_in is not None:        
        dMin_z = dMin_z_in
    else:
        dMin_z = np.nanmin(aData_reference-aData) #if it has negative value, change here  
    if (dMax_z <= dMin_z ):
        return
    aData_reference[np.where(aData_reference > dMax_z)] = dMax_z
    aData[np.where(aData > dMax_z)] = dMax_z #np.nan
    pv.set_plot_theme("document")
    
    surf = pv.StructuredGrid(aGrid_x, aGrid_y, aData_reference)
    wtd0 = pv.StructuredGrid(aGrid_x, aGrid_y, aData)
    wtd  = pv.StructuredGrid(aGrid_x, aGrid_y, aData_reference-aData)
    dRatio=1.0/100
    surf.points[:,2] *= dRatio
    wtd0.points[:,2] *= dRatio
    wtd.points[:,2] *= dRatio
    cmap0 = plt.get_cmap('terrain') 
    cmap1 = plt.get_cmap('rainbow')

    plotter = pv.Plotter()#off_screen=True)
    # Controlling the text properties
    sargs0 = dict(
        title_font_size=20,
        label_font_size=16,
        shadow=True,
        n_labels=5,
        italic=False,
        fmt="%.1f",
        font_family="arial",
        height=0.35, vertical=True,
        position_x=0.75, position_y=0.05,
        )
    sargs1 = dict(
        title_font_size=20,
        label_font_size=16,
        shadow=True,
        n_labels=5,
        italic=False,
        fmt="%.1f",
        font_family="arial",
        height=0.35, vertical=True,
        position_x=0.85, position_y=0.05,
        )

    #plotter.add_mesh(surf,scalars=surf.points[:,2]/dRatio,  stitle='Land surface (m)',  scalar_bar_args=sargs0, opacity=0.15)
    plotter.add_mesh(surf, show_edges=True, stitle='Land surface (m)',  color='white')
    
    #plotter.add_mesh(wtd, scalars=wtd0.points[:,2]/dRatio,  stitle='WTD (m)', cmap=cmap1, scalar_bar_args=sargs1)
    #below_color='blue', above_color='red')
   
    plotter.add_axes(xlabel='Longitude', ylabel='Latitude', zlabel='Elevation',line_width = 4)
    #plotter.enable_eye_dome_lighting()
    #plotter.enable_depth_peeling()
    
    #plotter.show_bounds( show_xaxis=True, show_yaxis=True, show_zaxis=True, \
    #    show_xlabels=True, show_ylabels=True, show_zlabels=False, italic=False, \
    #        bold=True, shadow=True, font_size=None, font_family=None, color=None, \
    #            xlabel='', ylabel='', zlabel='', \
    #                use_2d=False, grid=False, location='closest', ticks=None, \
    #                    all_edges=True, corner_factor=0.5, fmt=None, minor_ticks=False, padding=0.0)
    plotter.show_grid()

    #calcuale cpos based on data
    dMin_x = np.min(aGrid_x)
    dMax_x = np.max(aGrid_x)
    dRange_x = int(dMax_x - dMin_x)
    dMin_y = np.min(aGrid_y)
    dMax_y = np.max(aGrid_y)
    dRange_y = int(dMax_y - dMin_y)

    dMin_z = np.nanmin(aData_reference) * dRatio
    dMax_z = np.nanmax(aData_reference) * dRatio
    dRange_z = int(dMax_z - dMin_z)

    a = ( dMax_x + 1.2*dRange_x ,  dMin_y -  1.2 *dRange_y, dMax_z + 9 * dRange_z)
    b = ( dMin_x + 0.5*dRange_x ,  dMax_y -0.5 * dRange_y, dMin_z )
    #b = np.array( [dMin_x, dMax_y, dMin_z] ) #+ np.array( a)
    c = (0, 0, 1)

    cpos = [a,b,c]

    plotter.camera_position = cpos
    plotter.show( screenshot=sFilename_out )
    #plotter.screenshot(filename = sFilename_out, transparent_background=True)
   
    #
    cpos = plotter.camera_position
    print(cpos)
    print('finished')
    
    

   

   
