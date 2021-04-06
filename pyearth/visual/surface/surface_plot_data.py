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

sSystem_paths = os.environ['PATH'].split(os.pathsep)
sys.path.extend(sSystem_paths)
from pyearth.system.define_global_variables import *
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

def surface_plot_data_monthly(aData_all, \
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

    nstress = len(aTime)

    nslice = len(aData_all)
    if aColor_in is not None:
        aColor = aColor_in
    else:
        aColor = create_diverge_rgb_color_hex(nslice)

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
    ax.pbaspect = np.array([2.0, 10.0, 0.5])

    x1 = mdates.date2num( aTime )
    verts=[]
    ys = range(nslice)
    for iSlice in range(1, nslice+1):
        aData = aData_all[iSlice-1]
        nan_index = np.where(aData == missing_value)
        aData[nan_index] = np.nan
        good_index = np.where(  ~np.isnan(aData))
        #y1 = (aData[1])[0]
        y_top = aData
        xs = x1
        zs = aData
        verts.append(polygon_under_graph(xs, zs))
        pass

    poly = PolyCollection(verts, facecolors= aColor ,alpha=.6)
    ax.add_collection3d(poly, zs=ys, zdir='y')
    pYear = mdates.YearLocator(3)   # every year
    pMonth = mdates.MonthLocator()  # every month
    ax.axis('on')

    ax.xaxis._axinfo["grid"]['linewidth'] = 0.
    ax.yaxis._axinfo["grid"]['linewidth'] = 0.
    ax.zaxis._axinfo["grid"]['color'] = "grey"
    ax.zaxis._axinfo["grid"]['linestyle'] = "--"

    ax.xaxis.set_major_locator(pYear)
    #ax.xaxis.set_minor_locator(pMonth)
    sYear_format = mdates.DateFormatter('%Y')


    ax.xaxis.set_major_formatter(sYear_format)

    ax.set_xlabel(sLabel_x)
    ax.set_xlim3d(min(x1), max(x1))

    ax.set_ylim3d( nslice-1,0)
    ax.set_yticks( range(nslice) )
    aLabel_y =[]
    for i in aLabel_y_in:
        aLabel_y.append( i.title() )
        pass

    ax.set_yticklabels(aLabel_y,fontsize=13 )

    ax.set_ylabel(sLabel_y)

    if (iReverse_z==1):
        ax.set_zlim3d(dMax_z, dMin_z)
        pass
    else:
        ax.set_zlim3d(dMin_z, dMax_z)
        pass

    ax.set_zlabel(sLabel_z)




    plt.savefig(sFilename_out, bbox_inches='tight')

    plt.close('all')
    plt.clf()
    #print('finished plotting')
