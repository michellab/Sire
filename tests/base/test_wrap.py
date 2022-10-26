
import sire as sr

def _assert_array_equal( array0, array1 ):
    assert( len(array0) == len(array1) )

    for i in range(0, len(array0)):
        assert( array0[i] == array1[i] )

def test_wrap(verbose=False):
    from sire.base import wrap

    water = sr.load_test_files("water.pdb")[0]

    center = water.evaluate().center()

    dblarray = [ 1.0,2,3,4,5 ]
    intarray = [ 1,2,3,4,5 ]
    vecarray = [ sr.maths.Vector(1), sr.maths.Vector(2), sr.maths.Vector(3) ]
    strarray = [ "cat", "dog", "fish" ]
    x = sr.cas.Symbol("x")
    f = (x+5)**2
    mixarray = [ x, f, 5.3, "hello", sr.mol.Molecule(), [f, "cat", sr.vol.PeriodicBox()] ]

    water = water.edit().set_property("center", wrap(center)) \
                        .set_property("dblarray", wrap(dblarray)) \
                        .set_property("intarray", wrap(intarray)) \
                        .set_property("vecarray", wrap(vecarray)) \
                        .set_property("strarray", wrap(strarray)) \
                        .set_property("type", wrap("ligand")) \
                        .set_property("alpha", wrap(0.5)) \
                        .set_property("copies", wrap(1)) \
                        .set_property("mix", wrap(mixarray)).commit()

    assert water.property("center").value() == center
    _assert_array_equal( water.property("dblarray").value(), dblarray )
    _assert_array_equal( water.property("intarray").value(), intarray )
    _assert_array_equal( water.property("vecarray").value(), vecarray )
    _assert_array_equal( water.property("strarray").value(), strarray )
    assert water.property("type").value() == "ligand"
    assert water.property("alpha").value() == 0.5
    assert water.property("copies").value() == 1

    p = water.property("mix")

    assert sr.cas.Expression(x) == p[0].value() 
    assert f == p[1].value()
    assert p[2].value() == 5.3
    assert p[3].value() == "hello"
    assert p[4] == sr.mol.Molecule()
    assert p[5][0].value() == f
    assert p[5][1].value() == "cat"
    assert p[5][2] == sr.vol.PeriodicBox()

if __name__ == "__main__":
    test_wrap(True)
