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


import math
import time

from util import tinyr, node_gen, RectangleSetUniform, RectangleSetGaussian, RectangleSetCluster, RectangleSetMixed
import cProfile

def mean(values):
    return float(sum(values))/len(values)

def mean_and_deviation(items):
    m = mean(items)
    sqsum = 0
    for it in items:
        sqsum += (it-m)**2
    return (m, math.sqrt(sqsum/len(items)))

def create_datasets(entries=100000):
    return [('uniform', list(RectangleSetUniform(count=entries))),
            ('gaussian', list(RectangleSetGaussian(count=entries))),
            ('cluster', list(RectangleSetCluster(count=entries))),
            ('mixed', list(RectangleSetMixed(count=entries)))]
    

def timings():
    raise NotImplemented


if __name__ == '__main__':
    timings()


