
from Sire.MM import *
from Sire.Base import *

from nose.tools import assert_equal

def test_cast(verbose=False):
    c = StringProperty( VariantProperty("hello") )
    assert_equal( c, StringProperty("hello") )

    c = NumberProperty( VariantProperty(1.0) )
    assert_equal( c, NumberProperty(1.0) )

    c = NumberProperty( VariantProperty(5) )
    assert_equal( c, NumberProperty(5) )

    c = BooleanProperty( VariantProperty(False) )
    assert_equal( c, BooleanProperty(False) )

    c = BooleanProperty( VariantProperty(True) )
    assert_equal( c, BooleanProperty(True) )

def test_ff_cast(verbose=False):
    ff = InternalFF()
    ff.setProperty("combiningRules", VariantProperty("geometric"))


if __name__ == "__main__":
    test_cast(True)
    test_ff_cast(True)

