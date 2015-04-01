
from Sire.Base import *

from nose.tools import assert_equal

def test_refcount():
    mieow = StringProperty("mieow")
    p = Properties()
    p.setProperty("cat", mieow)

    #Deleting 'mieow' - should not crash!
    mieow = 0

    #accessing 'p' - should not crash!
    props = p.toString()

    mieow = p.property("cat")

    assert_equal(mieow, StringProperty("mieow"))

    mieow = 0

    #Deleting 'p' - should not crash!
    p = 0
