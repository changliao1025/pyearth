
import numpy as np

from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d
from pyearth.gis.geometry.calculate_angle_between_vectors_degrees import calculate_angle_between_vectors_degrees

def calculate_angle_betwen_vertex(dLongitude1_in, dLatitude1_in, 
                                dLongitude2_in, dLatitude2_in, 
                                dLongitude3_in, dLatitude3_in, 
                                iFlag_radian = None):
    #all in degree
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
    a3 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    b3 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    c3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)
    # Vectors in 3D space
    a3vec = a3 - b3
    c3vec = c3 - b3        
    angle3deg = calculate_angle_between_vectors_degrees(a3vec, c3vec)
    return  angle3deg



if __name__ == '__main__':

    dLongitude1_in=-148
    dLatitude1_in= 71
 
    dLongitude2_in=-148.1875
    dLatitude2_in=70.125
    
    dLongitude3_in=-147
    dLatitude3_in= 70.125
                          
    

    dAngle = calculate_angle_betwen_vertex(dLongitude1_in, dLatitude1_in, 
                                dLongitude2_in, dLatitude2_in, 
                                dLongitude3_in, dLatitude3_in 
                               )
    print(dAngle)