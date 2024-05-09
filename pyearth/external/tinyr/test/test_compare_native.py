

from testutils import make_tree_and_native, rand_rect, make_subrect, randrange_float


def test_search_superrect():
    rt, native, items = make_tree_and_native(howmuch=500, size=100, shift=800)
    
    for i in range(100):
        coords = rand_rect(size=200, shift=700)
        
        result = rt.search(coords)
        native_result = native.search(coords)
        assert set(result) == set(native_result), '%s; %s' % (str(result), str(native_result))

    

def test_interleaved():
    rt, native, items = make_tree_and_native(howmuch=200, size=70, shift=300, interleaved=True)
    
    for ident, coords in items:
        
        # exact search 
        result = rt.search(coords)
        native_result = native.search(coords)
        assert set(result) == set(native_result), '%s; %s' % (str(result), str(native_result))    
        
        # sub-rectangle search
        coords_subrect = make_subrect(coords, interleaved=True)
        result = rt.search(coords_subrect)
        native_result = native.search(coords_subrect)
        assert set(result) == set(native_result), '%s; %s' % (str(result), str(native_result))

        # point search
        point = [randrange_float(coords[0], coords[2]),
                 randrange_float(coords[1], coords[3]) ]
        
        result = rt.search_surrounding(point)
        native_result = native.search_surrounding(point)
        assert set(result) == set(native_result), '%s; %s' % (str(result), str(native_result))

if __name__ == '__main__':
    test_search_superrect()
    test_interleaved()
    print ('ok')
