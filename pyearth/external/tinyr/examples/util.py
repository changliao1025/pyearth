#
#    tinyr - a 2D-RTree implementation in Cython
#    Copyright (C) 2011  Matthias Simon
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


try:
    import tinyr
except ImportError as e:
    try:
        import pyximport
    except ImportError:
        raise e
    import sys, os
    sys.path.append(os.path.abspath(os.path.join('..', 'tinyr')))
    pyximport.install()
    import tinyr
    del pyximport

import random
import math
from collections import namedtuple

Node = namedtuple('Node', 'number coordinates')

###############################################################################
#--- Rectangle generators - ideas taken from [Bec90]
###############################################################################

class RectangleSet(object):
    '''Base for rectangle generators'''
    def __init__(self, *a, **kw):
        self.count = kw['count']
        self.interleaved = kw.get('interleaved', True)
        
    def __iter__(self):
        return self
    
    def next(self):
        raise NotImplementd

    def _distribute(self, a, b):
        # axial independent uniform distribution
        return (random.random() * (1-a), 
                random.random() * (1-b))

class RectangleSetUniform(RectangleSet):
    '''Rectangle generator with similar rectangles, randomly distributed'''
    def __init__(self, *a, **kw):
        RectangleSet.__init__(self, *a, **kw)
        self.u_area = 0.0001
        norm_var = 0.5
        self.deviation = self.u_area * norm_var
        
    def next(self):
        if self.count <= 0:
            raise StopIteration
        
        self.count -= 1
        
        area = -1
        while area <= 0:
            area = random.normalvariate(self.u_area, self.deviation)
            
        a = math.sqrt(area)
        
        # we don't want exact quadratic rectangles, so we re-shape
        a1 = -1
        while a1 <= 0: # we don't permit side with negative length
            a1 = a * random.normalvariate(1.0, 0.4)
        a = a1
        
        b = area / a

        x, y = self._distribute(a, b)

        x1 = x
        x2 = x + a
        y1 = y
        y2 = y + b
        
        if self.interleaved:
            return [x1, y1, x2, y2]
        return [x1, x2, y1, y2]

class RectangleSetGaussian(RectangleSetUniform):
    '''Rectangle generator with axial independent gaussian distribution'''
    def _distribute(self, a, b):
        while True:
            x = random.gauss(0.5, 0.2) - a/2
            y = random.gauss(0.5, 0.2) - b/2
            if x > 0 and y > 0:
                return x, y

class RectangleSetCluster(RectangleSetUniform):
    '''Rectangle generator with 16 clusters of rectangles'''
    def _distribute(self, a, b):
        cluser_row, cluser_col = divmod(random.randrange(0, 16), 4)
        
        # normal distribution with center on one of the 16 clusters
        x = random.normalvariate(0.5, 0.25) / 4 + cluser_col*0.25 - a/2
        y = random.normalvariate(0.5, 0.25) / 4 + cluser_row*0.25 - b/2
        
        return x , y
        
class RectangleSetMixed(RectangleSetUniform):
    '''Rectangle generator with 10% large and 90% small rectangles'''
    def __init__(self, *a, **kw):
        RectangleSetUniform.__init__(self, *a, **kw)
        
        self.state = 0
        self.u_area *= 10
        self.deviation = self.u_area * 0.5

        self.count *= 0.1
        kw['count'] *= 0.9
        
        self.small_rects = RectangleSetUniform(*a, **kw)
        
    def next(self):
        if self.state == 0:
            try:
                return self.small_rects.next()
            except StopIteration:
                self.state = 1
        
        return RectangleSetUniform.next(self)
        
###############################################################################

def rand_rect(size=100, shift=1000, interleaved=True):
    coords = [random.random() * size for i in range(4)]
    xshift = random.random() * shift
    yshift = random.random() * shift
    coords[0] += xshift
    coords[3] += yshift
    if interleaved:
        coords[1] += yshift
        coords[2] += xshift
        return coords
    
    else:
        coords[1] += xshift
        coords[2] += yshift
        return coords

def node_gen(count, size=100, shift=1000, interleaved=True, preordered_coordinates=False):
    if interleaved:
        idxs_off = 2
        y_start = 1
    else:
        idxs_off = 1
        y_start = 2
        
    for ident in range(count):
        coords = rand_rect(size, shift, interleaved)
        
        if preordered_coordinates:
            if coords[0] > coords[idxs_off]:
                coords[0], coords[idxs_off] = coords[idxs_off], coords[0]
            if coords[y_start] > coords[y_start+idxs_off]:
                coords[y_start], coords[y_start+idxs_off] = coords[y_start+idxs_off], coords[y_start]
                
        yield Node(ident, tuple(coords))


