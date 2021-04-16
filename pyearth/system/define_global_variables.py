#this module will be used to define all the global variable with consideration of cross platform
import os #retrieve existing system variables
import sys #used to add system path
import platform  #determine the platform this package is running at
from pathlib import Path #get the home directory
import getpass
sPlatform_os = platform.system()

sWorkspace_home = str(Path.home())


sUsername = getpass.getuser()


if sPlatform_os == 'Windows':  #windows
    slash = '\\'
    sMachine ='None'
       
    sWorkspace_scratch = 'C:'    
else:  #linux or unix
    slash = '/'    
   
    if (sPlatform_os == 'Linux'):
        sCluster = os.environ['SYSTEM_NAME']
        
        sWorkspace_scratch =  sWorkspace_home + slash + 'scratch'            
        
    else:
        if (sPlatform_os == 'Darwin'):
            sMachine ='mac'
            sWorkspace_scratch =  sWorkspace_home + slash + 'scratch'  
        else:
            pass

#system wide paths

#now we will start define major global variables
#data file type
sExtension_txt = '.txt'
sExtension_envi = '.dat'
sExtension_tiff = '.tif'
sExtension_header ='.hdr'
sExtension_netcdf = '.nc'
sExtension_shapefile = '.shp'

#graphics
sExtension_png = '.png'
sExtension_jpg = '.jpg'
sExtension_ps = '.ps'


#constant values
missing_value = -9999.0

nmonth = 12 #be careful with this one

#physical constants

iMonth_start = 1
iMonth_end = 12
mms2mmd = 24 * 3600.0
feet2meter = 0.3048
inch2mm = 25.4
cms2cmd = 24 * 3600
