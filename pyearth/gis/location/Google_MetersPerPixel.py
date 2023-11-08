import numpy as np
from pyearth.system.define_global_variables import *
def Google_MetersPerPixel( zoomLevel ):  
   
   # Return to the caller if there is an error.
   #On_Error, 2
   
   #; Need a zoom level?
   
   
   # Number of pixels in an image with a zoom level of 0.
   pixels_in_image = 256
   
   # The equitorial radius of the Earth assuming WGS-84 ellipsoid.
  
   
   # The number of meters per pixel.
   metersPerPixel = (2* np.pi * earth_radius) / pixels_in_image / np.power(2,zoomLevel)
   
   # Return the value.
   return metersPerPixel