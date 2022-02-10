
from math import cos, sin,  sqrt, pi
import numpy as np
def calculate_polygon_area(aLongitude_in, aLatitude_in,  iFlag_algorithm = 0, radius = 6378137.0):
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
    if npoint<3:
        print('More than 2 points are required!')
        return
   
    #TODO: take into account geodesy (i.e. convert latitude to authalic sphere, use radius of authalic sphere instead of mean radius of spherical earth)
    aLatitude_in = np.deg2rad(aLatitude_in)
    aLongitude_in = np.deg2rad(aLongitude_in)

    if iFlag_algorithm==0:
        # Line integral based on Green's Theorem, assumes spherical Earth
        

        #close polygon
        if aLatitude_in[0]!=aLatitude_in[-1]:
            aLatitude_in = np.append(aLatitude_in, aLatitude_in[0])
            aLongitude_in = np.append(aLongitude_in, aLongitude_in[0])

        # Get colatitude (a measure of surface distance as an angle)
        a = sin(aLatitude_in/2)**2 + cos(aLatitude_in)* sin(aLongitude_in/2)**2
        colat = 2*np.arctan2( sqrt(a), sqrt(1-a) )

        #azimuth of each point in segment from the arbitrary origin
        az = np.arctan2(cos(aLatitude_in) * sin(aLongitude_in), sin(aLatitude_in)) % (2*pi)

        # Calculate step sizes
        # daz = np.diff(az) % (2*pi)
        daz = np.diff(az)
        daz = (daz + pi) % (2 * pi) - pi

        # Determine average surface distance for each step
        deltas=np.diff(colat)/2
        colat=colat[0:-1]+deltas

        # Integral over azimuth is 1-cos(colatitudes)
        integrands = (1-cos(colat)) * daz

        # Integrate and save the answer as a fraction of the unit sphere.
        # Note that the sum of the integrands will include a factor of 4pi.
        area = abs(sum(integrands))/(4*pi) # Could be area of inside or outside

        area = min(area,1-area)
        if radius is not None: #return in units of radius
            return area * 4*pi*radius**2
        else: #return in ratio of sphere total area
            return area
    elif iFlag_algorithm==2:
        #L'Huilier Theorem, assumes spherical earth
        #see:
        # https://mathworld.wolfram.com/SphericalPolygon.html
        # https://web.archive.org/web/20160324191929/http://forum.worldwindcentral.com/showthread.php?20724-A-method-to-compute-the-area-of-a-spherical-polygon
        # https://github.com/spacetelescope/spherical_geometry/blob/master/spherical_geometry/polygon.py
        # https://github.com/tylerjereddy/spherical-SA-docker-demo/blob/master/docker_build/demonstration.py
        #TODO
        pass
    elif iFlag_algorithm==3:
        #https://trs.jpl.nasa.gov/handle/2014/41271
        #TODO
        pass