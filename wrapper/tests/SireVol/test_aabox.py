
from Sire.Vol import *
from Sire.Maths import *

from nose.tools import assert_almost_equal

def test_within(verbose=False):
    aabox0 = AABox(Vector(0), Vector(3,4,5))

    aabox1 = AABox(Vector(10,0,0), Vector(2,2,2))

    assert( aabox0.withinDistance(5.1, aabox1) )

def test_space(verbose=False):
    space = Cartesian()

    aabox0 = AABox(Vector(0), Vector(3,4,5))
    aabox1 = AABox(Vector(10,0,0), Vector(2,2,2))

    rand = RanGenerator()

    for i in range(0,100):
        delta = Vector( rand.rand(-5,5), rand.rand(-5,5), rand.rand(-5,5) )

        aabox0.translate(delta)
        aabox1.translate(delta) 

        mindist0 = space.minimumDistance(aabox0, aabox1)
        mindist1 = space.minimumDistance(aabox1, aabox0)

        if verbose:
            print("%s %s %s %s" % (aabox0, aabox1, mindist0, mindist1))

        assert_almost_equal(mindist0, 5)
        assert_almost_equal(mindist1, 5)


if __name__ == "__main__":
    test_within(True)
    test_space(True)
