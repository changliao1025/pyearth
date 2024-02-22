

import numpy as np
from pyearth.gis.location.convert_longitude_latitude_to_sphere_3d import convert_longitude_latitude_to_sphere_3d

def calculate_angle_between_vertex_normal(dLongitude1_in, dLatitude1_in, 
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
    a3 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    b3 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    c3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)
    a3vec = a3 - b3
    c3vec = c3 - b3 
    dot = np.dot(a3vec, c3vec)
    g = np.cross(a3vec, c3vec)
    det = np.dot(b3, g)
    angle = np.arctan2(det, dot)
    f = np.degrees(angle) 
    if f < 0:
        f = 360 + f
    
    return f


if __name__ == '__main__':

    dLongitude1_in=-149
    dLatitude1_in= 71
 
    dLongitude2_in=-148.1875
    dLatitude2_in=70.125
    
    dLongitude3_in=-147
    dLatitude3_in= 70.125
                          
    

    dAngle = calculate_angle_between_vertex_normal(dLongitude1_in, dLatitude1_in, 
                                dLongitude2_in, dLatitude2_in, 
                                dLongitude3_in, dLatitude3_in 
                               )
    print(dAngle)


    #calculate left or right of line segment

    #a
    dLongitude1_in=-149
    dLatitude1_in=70

    #b
    dLongitude2_in=-148
    dLatitude2_in= 70   
    
    #c
    dLongitude3_in=-147
    dLatitude3_in= 70

    #d
    dLongitude4_in=-149
    dLatitude4_in= 72

    dLongitude1_radian_in, dLatitude1_radian_in = np.radians(np.array((dLongitude1_in, dLatitude1_in) ))
    dLongitude2_radian_in, dLatitude2_radian_in = np.radians(np.array((dLongitude2_in, dLatitude2_in) )) #this is the middle one
    dLongitude3_radian_in, dLatitude3_radian_in = np.radians(np.array((dLongitude3_in, dLatitude3_in) ))
    dLongitude4_radian_in, dLatitude4_radian_in = np.radians(np.array((dLongitude4_in, dLatitude4_in) ))
    a3 = convert_longitude_latitude_to_sphere_3d(dLongitude1_radian_in, dLatitude1_radian_in)
    b3 = convert_longitude_latitude_to_sphere_3d(dLongitude2_radian_in, dLatitude2_radian_in)
    c3 = convert_longitude_latitude_to_sphere_3d(dLongitude3_radian_in, dLatitude3_radian_in)  
    d3 = convert_longitude_latitude_to_sphere_3d(dLongitude4_radian_in, dLatitude4_radian_in)

    a3vec = b3 - a3
    c3vec = c3 - b3 
    d3vec = b3 - d3


    dot = np.dot(c3vec, a3vec)
    g = np.cross(c3vec, a3vec)
    det = np.dot(b3, g)
    angle1 = np.arctan2(det, dot)
    f1 = np.degrees(angle1) 
    if f1 < 0:
        f1 = 360 + f1
    

    dot = np.dot(d3vec, a3vec)
    g = np.cross(d3vec, a3vec)
    det = np.dot(b3, g)
    angle2 = np.arctan2(det, dot)
    f2 = np.degrees(angle2) 
    if f2 < 0:
        f2 = 360 + f2
                          
    
    print(f1)
    print(f2)

    delta = f1 - f2 
    print(delta)
    