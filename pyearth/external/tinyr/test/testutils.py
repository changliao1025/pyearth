
import os
import sys
import random

tinyr_dir = os.path.dirname(os.path.dirname(__file__))

p = os.path.abspath(os.path.join(tinyr_dir, 'examples'))
sys.path.append(p)
p = os.path.abspath(os.path.join(tinyr_dir, 'tinyr'))
sys.path.append(p)

from util import tinyr, node_gen, rand_rect

def randrange_float(start, stop):
    return random.random() * (stop-start) + start

def make_tree(howmuch=100, size=100, shift=800, interleaved=False, preordered_coordinates=False):
    rt = tinyr.RTree(interleaved=interleaved, max_cap=5, min_cap=2)
    nodegen = node_gen(howmuch, size=size, shift=shift, interleaved=interleaved, preordered_coordinates=preordered_coordinates)
    
    items = []
    for item in nodegen:
        items.append(item)
        
        ident, coords = item
        rt.insert(ident, coords)

        rt.valid()

    return rt, items

def make_tree_and_native(howmuch=100, size=100, shift=800, interleaved=False):
    rt = tinyr.RTree(interleaved=interleaved, max_cap=5, min_cap=2)
    native = FakeTree(interleaved=interleaved)
    nodegen = node_gen(howmuch, size=size, shift=shift, interleaved=interleaved)
    
    items = []
    for item in nodegen:
        items.append(item)
        
        ident, coords = item
        rt.insert(ident, coords)
        native.insert(ident, coords)
        rt.valid()

    return rt, native, items


def make_subrect(coords, interleaved=False):
    # random shrink factor
    if interleaved:
        maxshrink_x = (coords[2] - coords[0]) / 2.0 - 1
        maxshrink_y = (coords[3] - coords[1]) / 2.0 - 1
    else:
        maxshrink_x = (coords[1] - coords[0]) / 2.0 - 1
        maxshrink_y = (coords[3] - coords[2]) / 2.0 - 1
    
    coords = list(coords)
    coords[0] += randrange_float(0, maxshrink_x)
    coords[3] -= randrange_float(0, maxshrink_y)
    
    if interleaved:
        coords[1] += randrange_float(0, maxshrink_y)
        coords[2] -= randrange_float(0, maxshrink_x)
    else:
        coords[2] += randrange_float(0, maxshrink_y)
        coords[1] -= randrange_float(0, maxshrink_x)
    
    return coords

class FakeTree(object):
    '''An index structure with the API of r-tree - to compare results'''
    def __init__(self, interleaved=False):
        self.items = []
        self.interleaved = interleaved
    
    def _overlaps(self, a, b):
        return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]

    def _make_order(self, coords):
        coords = list(coords)
        
        if not self.interleaved:
            coords[1], coords[2] =  coords[2], coords[1]
            
        if coords[0] > coords[2]:
            coords[0], coords[2] = coords[2], coords[0]
        if coords[1] > coords[3]:
            coords[1], coords[3] = coords[3], coords[1]
            
        return tuple(coords)
    
    def insert(self, ident, coords):
        coords = self._make_order(coords)
        assert coords[0] <= coords[2] and  coords[1] <= coords[3], coords
        self.items.append((ident, coords))
    
    def search(self, coords):
        coords = self._make_order(coords)
        return [ ident for ident, rect in self.items if self._overlaps(coords, rect) ]
    
    def search_surrounding(self, point):
        coords = point*2
        coords = self._make_order(coords)
        return [ ident for ident, rect in self.items if self._overlaps(coords, rect) ]
        
    def remove(self, coords):
        coords = self._make_order(coords)
            
        for item in self.items:
            ident, rect = item
            if self._overlaps(coords, rect):
                self.items.remove(item)
                break
        else:
            raise KeyError


