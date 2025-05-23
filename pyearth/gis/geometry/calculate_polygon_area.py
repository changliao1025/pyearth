#it is also recommended to use the nvector API to calculate the area of a polygon on the ellipsoid
import math
import numpy as np
from pyearth.system.define_global_variables import *

from pyearth.gis.geometry.calculate_spherical_triangle_area import calculate_spherical_triangle_area
from pyearth.gis.geometry.calculate_distance_based_on_longitude_latitude import calculate_distance_based_on_longitude_latitude


def haversine(x):
    """
    Haversine function: hav(x) = (1 - cos(x)) / 2
    :param x: Angle in radians
    :return: Returns the value of the Haversine function
    """
    return (1.0 - math.cos(x)) / 2.0

def calculate_polygon_area(aLongitude_in,
                           aLatitude_in,
                           iFlag_algorithm=1,
                           iFlag_radian=None,
                           dRadius_in=None,
                           dLine_threshold = None):
    """
    Computes area of spherical polygon, assuming spherical Earth.
    Returns result in ratio of the sphere's area if the radius is specified. Otherwise, in the units of provided radius.
    lats and lons are in degrees.

    Args:
        aLongitude_in (list): The longitude of list of points
        aLatitude_in (list): The latitude of list of points
        iFlag_algorithm (int, optional): Which algorithm to use. Defaults to 0.
        radius (float, optional): The radius of Earth in meter. Defaults to 6378137.0.

    Returns:
        float: The area of the polygon
    """

    npoint = len(aLongitude_in)
    if npoint < 3:
        print('More than 2 points are required!')
        return

    # TODO: take into account geodesy (i.e. convert latitude to authalic sphere, use radius of authalic sphere instead of mean radius of spherical earth)
    # close polygon
    if aLatitude_in[0] != aLatitude_in[-1] or aLongitude_in[0] != aLongitude_in[-1]:
        aLatitude_in = np.append(aLatitude_in, aLatitude_in[0])
        aLongitude_in = np.append(aLongitude_in, aLongitude_in[0])

    npoint = len(aLongitude_in)

    if iFlag_radian is None:  # degree_based
        aLongitude_radian_in = np.deg2rad(aLongitude_in)
        aLatitude_radian_in = np.deg2rad(aLatitude_in)

    else:
        aLongitude_radian_in = aLongitude_in
        aLatitude_radian_in = aLatitude_in
        pass



    #check whether the polygon is close to line, narrow and thin:
    # Calculate the lengths of the sides of the polygon
    aLength = np.zeros(npoint-1)
    for i in range(npoint-1):
        dLength = calculate_distance_based_on_longitude_latitude(aLongitude_in[i], aLatitude_in[i], aLongitude_in[i+1], aLatitude_in[i+1])
        aLength[i] = dLength
        pass

    dLength_max = np.max(aLength)
    dlength_rest = np.sum(aLength) - dLength_max
    # Check if the polygon is close to a line
    if dLine_threshold is not None:
        if (dLength_max / dlength_rest) > (1 - dLine_threshold):
            print( "The polygon is close to a line" )
            area = 0.0
            return area
        else:
            pass

    if iFlag_algorithm == 0:
        # Line integral based on Green's Theorem, assumes spherical Earth

        # Get colatitude (a measure of surface distance as an angle)
        a = np.sin(aLatitude_radian_in/2)**2 + np.cos(aLatitude_radian_in) * \
            np.sin(aLongitude_radian_in/2)**2
        colat = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))

        # azimuth of each point in segment from the arbitrary origin
        az = np.arctan2(np.cos(aLatitude_radian_in) * np.sin(aLongitude_in),
                        np.sin(aLatitude_radian_in)) % (2*np.pi)

        # Calculate step sizes
        # daz = np.diff(az) % (2*pi)
        daz = np.diff(az)
        daz = (daz + np.pi) % (2 * np.pi) - np.pi

        # Determine average surface distance for each step
        deltas = np.diff(colat)/2
        colat = colat[0:-1]+deltas

        # Integral over azimuth is 1-cos(colatitudes)
        integrands = (1-np.cos(colat)) * daz

        # Integrate and save the answer as a fraction of the unit sphere.
        # Note that the sum of the integrands will include a factor of 4pi.
        area = abs(sum(integrands))  # Could be area of inside or outside

        area = min(area, 1-area)

    elif iFlag_algorithm == 1:
        # L'Huilier Theorem, assumes spherical earth
        # see:
        # https://mathworld.wolfram.com/SphericalPolygon.html
        # https://web.archive.org/web/20160324191929/http://forum.worldwindcentral.com/showthread.php?20724-A-method-to-compute-the-area-of-a-spherical-polygon
        # https://github.com/spacetelescope/spherical_geometry/blob/master/spherical_geometry/polygon.py
        # https://github.com/tylerjereddy/spherical-SA-docker-demo/blob/master/docker_build/demonstration.py
        #

        dLongtitude_root = aLongitude_radian_in[0]
        dLatitude_root = aLatitude_radian_in[0]
        dLongtitude_b = aLongitude_radian_in[1]
        dLatitude_b = aLatitude_radian_in[1]
        nTriangle = npoint - 2
        aArea = np.zeros(nTriangle)
        for i in np.arange(1, nTriangle+1, 1):
            # define a triangle
            dLongtitude_a = dLongtitude_b
            dLatitude_a = dLatitude_b
            dLongtitude_b = aLongitude_radian_in[i+1]
            dLatitude_b = aLatitude_radian_in[i+1]
            # calculate the area of the triangle
            aLongitude_temp = [dLongtitude_root, dLongtitude_a, dLongtitude_b]
            aLatitude_temp = [dLatitude_root, dLatitude_a, dLatitude_b]
            dArea_triangle = calculate_spherical_triangle_area(aLongitude_temp,
                                                               aLatitude_temp,
                                                               iFlag_radian=1)

            aArea[i-1] = dArea_triangle
            pass

        area = np.sum(aArea)
        pass

    elif iFlag_algorithm == 2:
        # https://trs.jpl.nasa.gov/handle/2014/41271
        # TODO
        if dRadius_in is not None:
            dRadius = dRadius_in
        else:
            dRadius = earth_radius
        dArea_m = spherical_polygon_area(aLatitude_radian_in, aLongitude_radian_in, dRadius)
        return float(dArea_m)
        pass

    if iFlag_radian is not None:
        return float(area)
    else:
        #6371229.0 or something else
        if dRadius_in is not None:
            dArea_m = area * dRadius_in**2
        else:
            dArea_m = area * earth_radius**2

        return float(dArea_m)


def spherical_polygon_area(lat, lon, r):
    """
    Compute the Area of a Spherical Polygon
    :param lat: List of latitudes of all vertices (in radians)
    :param lon: List of longitudes of all vertices (in radians)
    :param r: Spherical radius
    :return: Returns the area of a spherical polygon
    """
    lam1 = lam2 = beta1 = beta2 = cosB1 = cosB2 = 0
    hav = 0
    sum = 0

    for j in range(len(lat)):
        k = j + 1
        if j == 0:
            lam1 = lon[j]
            beta1 = lat[j]
            lam2 = lon[j + 1]
            beta2 = lat[j + 1]
            cosB1 = math.cos(beta1)
            cosB2 = math.cos(beta2)
        else:
            k = (j + 1) % len(lat)
            lam1 = lam2
            beta1 = beta2
            lam2 = lon[k]
            beta2 = lat[k]
            cosB1 = cosB2
            cosB2 = math.cos(beta2)

        if lam1 != lam2:
            hav = haversine(beta2 - beta1) + cosB1 * cosB2 * haversine(lam2 - lam1)
            a = 2 * math.asin(math.sqrt(hav))
            b = math.pi / 2 - beta2
            c = math.pi / 2 - beta1
            s = 0.5 * (a + b + c)
            t = math.tan(s / 2) * math.tan((s - a) / 2) * math.tan((s - b) / 2) * math.tan((s - c) / 2)

            excess = abs(4 * math.atan(math.sqrt(abs(t))))

            if lam2 < lam1:
                excess = -excess

            sum += excess

    return abs(sum) * r * r

if __name__ == '__main__':
    # test the polygon area calculation
    #aLongitude_in = [-148, -148, -148.5, -148.5]
    #aLatitude_in = [68, 68.5, 68.5, 68]
#
    #dArea = calculate_polygon_area(aLongitude_in, aLatitude_in, iFlag_algorithm=0)
    #print(dArea)
#
    #dArea1 = calculate_polygon_area(aLongitude_in, aLatitude_in)
    #print(dArea1)
#
    #print(dArea/dArea1)

    xc = -148.1875
    yc = 70.125
    dResolution = 1.0/16

    aLongitude_in = [xc-dResolution/2, xc+dResolution/2, xc+dResolution/2, xc-dResolution/2]
    aLatitude_in = [yc-dResolution/2, yc-dResolution/2, yc+dResolution/2, yc+dResolution/2]
    dArea1 = calculate_polygon_area(aLongitude_in, aLatitude_in)
    print(dArea1)
