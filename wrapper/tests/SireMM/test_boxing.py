
from Sire.MM import *
from Sire.Maths import *
from Sire.Units import *

def test_boxing():
    length = 5 * angstrom

    for ix in range(-200,200,1):
        x = 0.1 * ix

        for iy in range(-200,200,13):
            y = 0.1 * iy

            for iz in range(-200,200,23):
                z = 0.1 * iz

                v = Vector(x,y,z)

                cljindex = CLJBoxIndex.createWithBoxLength( v, length )

                aabox = cljindex.box(length)

                if not aabox.contains(v):
                    print("%s contains %s? %s %s" % (aabox.toString(), \
                              v.toString(), cljindex.toString(), aabox.contains(v)) )

                assert( aabox.contains(v) )
                

if __name__ == "__main__":
    test_boxing()
