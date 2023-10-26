import os
import numpy as np

import cartopy.crs as ccrs
import cartopy.mpl.ticker as ticker
import matplotlib as mpl
#from shapely.wkt import loads
from osgeo import  osr, gdal, ogr

import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.gis.gdal.gdal_functions import get_geometry_coords
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

pProjection = ccrs.PlateCarree() #for latlon data only

class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format

def map_netcdf_file(
                    sFilename_netcdf_in, 
                             
                             sFilename_out,    
                             
                               iFlag_unstructured_in =None,
                             sVariable_in = None,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFont_size_in = None,
                            sColormap_in = None,
                            sTitle_in = None, 
                            iDPI_in = None,
                            iSize_x_in = None, 
                            iSize_y_in = None, 
                            dMissing_value_in=None,
                            dData_max_in = None, 
                            dData_min_in = None,
                            sExtend_in =None,
                            sFont_in = None,
                            sUnit_in=None,
                            aLegend_in = None,
                            aExtent_in = None,
                            pProjection_map_in=None,
                               aBoundary_in=None):
    """
    Extract data from a NetCDF file based on a bounding box.

    Args:
        sFilename_netcdf_in (str): Path to the input NetCDF file.
        sFilename_netcdf_out (str): Path to the output NetCDF file.
        bbox (tuple): Bounding box as a tuple in the format (aLongitude_box_min, aLatitude_box_min, aLongitude_box_max, aLatitude_box_max).
    """
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iSize_x_in is not None:
        iSize_x = iSize_x_in
    else:
        iSize_x = 8

    if iSize_y_in is not None:
        iSize_y = iSize_y_in
    else:
        iSize_y = 8
    
    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

    if dData_min_in is not None:    
        iFlag_data_min = 1
        dData_min = dData_min_in
    else:
        iFlag_data_min = 0
        pass

    if dData_max_in is not None:    
        iFlag_data_max = 1
        dData_max = dData_max_in
    else:        
        iFlag_data_max = 0
        pass


    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'rainbow'

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title =1
    else:
        iFlag_title=0
        sTitle =  ''

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend =  'max'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit =  ''
    
    if sFont_in is not None:
        sFont = sFont_in
    else:    
        sFont = "Times New Roman"

    plt.rcParams["font.family"] = sFont
    
    if sVariable_in is not None:
        sVariable = sVariable_in
    else:
        sVariable =  'id'

    cmap = cm.get_cmap(sColormap)

    fig = plt.figure( dpi = iDPI  )
    
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    pSrs = osr.SpatialReference()  
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
   

    try:
        ds = nc.Dataset(sFilename_netcdf_in, 'r')
        aLatitude = ds.variables['lat'][:]
        aLongitude = ds.variables['lon'][:]

        aLongitude_box_min, aLatitude_box_min, aLongitude_box_max, aLatitude_box_max = bbox

        if iFlag_unstructured_in is not None:
            #unstructured grid
            pass
        else:   
            #regular grid
            dLongitude_min = np.min(aLongitude)
            dLongitude_max = np.max(aLongitude)
            dLatitude_min = np.min(aLatitude)
            dLatitude_max = np.max(aLatitude)
            ncolum = len(aLongitude)
            nrow = len(aLatitude)
            dResolution_x = (dLongitude_max - dLongitude_min) / (ncolum - 1)
            dResolution_y = (dLatitude_max - dLatitude_min) / (nrow - 1)

            #the bounding box should be applied to cell center

            for i in range(0, ncolum):
                for j in range(0, nrow):
                    dLongitude_cell_center = dLongitude_min + i * dResolution_x
                    dLatitude_cell_center = dLatitude_min + j * dResolution_y
                    if (dLongitude_cell_center >= aLongitude_box_min) \
                    and (dLongitude_cell_center <= aLongitude_box_max) \
                    and (dLatitude_cell_center >= aLatitude_box_min)\
                    and (dLatitude_cell_center <= aLatitude_box_max):
                        
                        pass
                    else:
                        pass
                    
        

        print("Extraction successful!")
    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage:

if __name__ == "__main__":
    sFilename_netcdf_in_path = "path_to_input_netcdf_file.nc"
    sFilename_netcdf_out_path = "path_to_output_netcdf_file.nc"
    bbox = (-90, 30, -70, 40)  # Example bounding box

    map_netcdf_file(sFilename_netcdf_in_path, sFilename_netcdf_out_path, bbox)
