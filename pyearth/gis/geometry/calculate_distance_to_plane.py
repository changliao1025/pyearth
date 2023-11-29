import numpy as np
from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d

#https://stackoverflow.com/questions/8204998/how-to-check-if-a-pointlonc-latc-lie-on-a-great-circle-running-from-lona-lata

def calculate_distance_to_plane(dLongitude1_in, dLatitude1_in, 
                                dLongitude2_in, dLatitude2_in, 
                                dLongitude3_in, dLatitude3_in, 
                                iFlag_radian = None):
    
    if iFlag_radian is None:
        dLongitude1_radian_in, dLatitude1_radian_in = np.radians(np.array((dLongitude1_in, dLatitude1_in) ))
        dLongitude2_radian_in, dLatitude2_radian_in = np.radians(np.array((dLongitude2_in, dLatitude2_in) )) #this is the middle one
        dLongitude3_radian_in, dLatitude3_radian_in = np.radians(np.array((dLongitude3_in, dLatitude3_in) ))
        pass
    else:         
        dLongitude1_radian_in, dLatitude1_radian_in = dLongitude1_in , dLatitude1_in  
        dLongitude2_radian_in, dLatitude2_radian_in = dLongitude2_in , dLatitude2_in #this is the middle one
        dLongitude3_radian_in, dLatitude3_radian_in = dLongitude3_in , dLatitude3_in
        pass
        
    # The points in 3D space
    x1,y1,z1 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    x2,y2,z2 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    x3,y3,z3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)
    #The formula is x+b*y+c*z=0 
    
    c = (-x1*y3 + x3* y1)/( z1*y3 - z3*y1 )
    b = (  -x1*z3 + x3 * z1 ) / (y1 * z3 - y3*z1)
    distance = abs(  x2 + b * y2 + c * z2 )
    return distance