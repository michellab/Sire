
from Sire.IO import *
from Sire.Mol import *
from Sire.Maths import *
from Sire.Base import *

from nose.tools import assert_equal

def assert_array_equal( array0, array1 ):
    assert( len(array0) == len(array1) )

    for i in range(0, len(array0)):
        assert( array0[i] == array1[i] )

def test_wrap():
    water = PDB().readMolecule("../io/water.pdb")

    center = water.evaluate().center()

    dblarray = [ 1.0,2,3,4,5 ]
    intarray = [ 1,2,3,4,5 ]
    vecarray = [ Vector(1), Vector(2), Vector(3) ]
    strarray = [ "cat", "dog", "fish" ]

    water = water.edit().setProperty("center", wrap(center)) \
                        .setProperty("dblarray", wrap(dblarray)) \
                        .setProperty("intarray", wrap(intarray)) \
                        .setProperty("vecarray", wrap(vecarray)) \
                        .setProperty("strarray", wrap(strarray)) \
                        .setProperty("type", wrap("ligand")) \
                        .setProperty("alpha", wrap(0.5)) \
                        .setProperty("copies", wrap(1)).commit()

    assert_equal( water.property("center"), center )
    assert_array_equal( water.property("dblarray"), dblarray )
    assert_array_equal( water.property("intarray"), intarray )
    assert_array_equal( water.property("vecarray"), vecarray )
    assert_array_equal( water.property("strarray"), strarray )
    assert_equal( water.property("type"), "ligand" )
    assert_equal( water.property("alpha"), 0.5 )
    assert_equal( water.property("copies"), 1 )
