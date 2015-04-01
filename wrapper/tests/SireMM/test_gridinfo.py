
from Sire.MM import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Units import *


def test_index(verbose = False):
    info = GridInfo( AABox( Vector(0), Vector(2) ), 1.0*angstrom )

    if verbose:
        print("Test grid: %s" % info)

    for i in range(0,info.nPoints()):
        idx = info[i]
        if verbose:
            print("%s: %s %s" % (i,idx,info.point(i)))
    
        assert( i == info[idx] )

    for i in range(0, info.nPoints()):
        if verbose:
            print("%s contains %s?" % (info.box(i), info.point(i)))

        assert( info.box(i).contains(info.point(i)) )

        if i < info.nPoints() - 1:
            assert( not info.box(i).contains(info.point(i+1)) )

def test_setup(verbose = False):
    info = GridInfo( AABox( Vector(0), Vector(2) ), 1.0*angstrom )

    if verbose:
        print("Test grid: %s" % info)

    assert( info.dimensions() == AABox( Vector(0), Vector(2) ) )
    assert( info.spacing() == 1*angstrom )
    assert( info.nPoints() == 125 )
    assert( info.dimX() == 5 )
    assert( info.dimY() == 5 )
    assert( info.dimZ() == 5 )
    assert( info.isEmpty() != True )

if __name__ == "__main__":
    test_setup(True)
    test_index(True)

