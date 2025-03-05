import numpy as np
from pyearth.gis.geometry.calculate_distance_to_line import calculate_distance_to_line

def douglas_peucker_geodetic(aCoords_gcs, tolerance):
    """
    Simplifies a polyline using the Douglas-Peucker algorithm.

    Args:
        aPoint (list of tuples): The input polyline as a list of (x, y) tuples.
        tolerance (float): The tolerance for simplification.

    Returns:
        list of tuples: The simplified polyline.
    """
    aPoint = list()
    for aCoord in aCoords_gcs:
        aPoint.append( (aCoord[0], aCoord[1]) )

    is_polygon = aPoint[0] == aPoint[-1]
    if is_polygon:
        aPoint = aPoint[:-1]
    # Find the point with the maximum distance from the line segment
    def find_furthest_point(points, start, end):
        max_distance = 0
        index = 0
        for i in range(start + 1, end):
            distance = point_line_distance(points[i], points[start], points[end])
            if distance > max_distance:
                max_distance = distance
                index = i
        return index, max_distance

    # Calculate the perpendicular distance from a point to a line segment
    def point_line_distance(point, start, end):
        #get the longitude and latitude from each point
        #if start == end:
        #    return np.linalg.norm(np.array(point) - np.array(start))
        #else:
        #    n = np.abs((end[1] - start[1]) * point[0] - (end[0] - start[0]) * point[1] + end[0] * start[1] - end[1] * start[0])
        #    d = np.linalg.norm(np.array(end) - np.array(start))
        #    return n / d
        dLongitude0 = point[0]
        dLatitude0 = point[1]
        dLongitude1 = start[0]
        dLatitude1 = start[1]
        dLongitude2 = end[0]
        dLatitude2 = end[1]
        #dLongitude0 dLatitude0 is the middle point

        distance = calculate_distance_to_line(
                dLongitude1, dLatitude1, dLongitude0, dLatitude0, dLongitude2, dLatitude2)
        #print(distance)
        return distance

    # Recursive function to simplify the polyline or polygon
    def simplify(points, start, end, tolerance, simplified):
        index, max_distance = find_furthest_point(points, start, end)
        if max_distance > tolerance:
            simplify(points, start, index, tolerance, simplified)
            simplify(points, index, end, tolerance, simplified)
        else:
            if start not in simplified:
                simplified.append(start)
            if end not in simplified:
                simplified.append(end)
    simplified = []
    simplify(aPoint, 0, len(aPoint) - 1, tolerance, simplified)
    simplified.sort()
    simplified_points = [aPoint[i] for i in simplified]
    # Close the polygon if it was originally closed
    if is_polygon:
        simplified_points.append(simplified_points[0])

    return simplified_points


# Example usage

if __name__ == '__main__':

    #create a polygon using a list of longitude and latitude location



    pass