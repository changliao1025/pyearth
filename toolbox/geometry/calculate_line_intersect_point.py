import shapely
from shapely.geometry import LineString, Point

#find the intersection of two lines

def calculate_line_intersect_point(pPointA, pPointB, pPointC, pPointD):
    pLine1 = LineString([pPointA, pPointB])
    pLine2 = LineString([pPointC, pPointD])

    pIntersect = pLine1.intersection(pLine2)
    pPoint_intersect = pIntersect.x, pIntersect.y
    
    print(pPoint_intersect)
    return pPoint_intersect 