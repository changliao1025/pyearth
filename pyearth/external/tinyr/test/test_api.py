


import StringIO

from testutils import make_tree, make_subrect, randrange_float, tinyr


def test_insert():
    make_tree()

def test_valid_empty():
    r = tinyr.RTree(min_cap=3, max_cap=7)
    assert r.valid()
    
    assert r.min_cap == 3
    assert r.max_cap == 7
    

def test_search_exact_match():
    rt, items = make_tree()
    for ident, coords in items:
        result = rt.search(coords)
        assert ident in result, (ident, result, coords)
    
    assert rt.valid()
    
def test_search_subrect():
    rt, items = make_tree()
    
    for ident, coords in items:
        coords = make_subrect(coords)
        result = rt.search(coords)
        assert ident in result, (ident, result, coords)

    assert rt.valid()

def test_search_point():
    rt, items = make_tree()
    
    for ident, coords in items:
        # random point in rect defined by coords
        point = [randrange_float(coords[0], coords[1]),
                 randrange_float(coords[2], coords[3])]
        
        result = rt.search_surrounding(point)
        assert ident in result, (ident, result, coords)
        
    assert rt.valid()


def test_remove():
    rt, items = make_tree()
    
    for item in items:
        rt.remove(*item)
        rt.valid()
    
    assert len(rt) == 0


def test_copy_empty():
    r1 = tinyr.RTree()
    r2 = r1.copy()
    
    assert len(r1) == len(r2)

def test_copy_full():
    rt, items = make_tree()
    rt.valid()
    
    rtnew = rt.copy()
    
    assert len(rtnew) == len(rt)
    # TODO: more equality checks
    assert rtnew.valid()

def test_interleaved():
    rt, items = make_tree(howmuch=200, size=70, shift=300, interleaved=True)
    
    for ident, coords in items:
        
        # exact search 
        result = rt.search(coords)
        assert ident in result, (ident, result, coords)
        
        # sub-rectangle search
        coords_subrect = make_subrect(coords, interleaved=True)
        result = rt.search(coords_subrect)
        assert ident in result, (ident, result, coords)

        # point search
        point = [randrange_float(coords[0], coords[2]),
                 randrange_float(coords[1], coords[3]) ]
        
        result = rt.search_surrounding(point)
        assert ident in result, (ident, result, coords)

def test_iter_kvi():
    
    rt, items = make_tree(preordered_coordinates=True)
    
    keys = [ n.coordinates for n in items ]
    values = [ n.number for n in items ]
    items = list(items)
    
    for k in rt.iterkeys():
        # if we wouldn't feed RTree with pre-ordered coordinates 
        # (x1<=x2 and y1<=y2), this would fail!
        keys.remove(k)
    
    assert len(keys) == 0, len(keys)
    
    for v in rt.itervalues():
        values.remove(v)
    
    assert len(values) == 0, len(values)

    for k, v in rt.iteritems():
        # again, this would fail with unordered coordinates
        i = (v, k)
        items.remove(i)
    
    assert len(items) == 0, len(items)
    

def test_infoobj():
    rt, items = make_tree(howmuch=100)
    
    for rt in (rt, tinyr.RTree()):
    
        io = rt.get_info()
        
        assert io.common_boundary
        assert io.width
        assert io.height
        
        assert len(list(io.iter_rectangles())) > 0
        
        io.to_dot(StringIO.StringIO())
        
