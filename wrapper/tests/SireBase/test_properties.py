
from Sire.Base import *

from nose.tools import assert_equal

def test_set_property():
    p = Properties()

    p.setProperty("author", wrap("Christopher"))

    assert_equal( p.property("author"), wrap("Christopher") )

def test_set_metadata():
    p = Properties()

    p.setMetadata("about", wrap([1,2,3,4]))

    assert_equal( p.metadata("about"), wrap([1,2,3,4]) )

