
import pytest

from sire.units import fahrenheit, celsius, farad, meter, convert
from sire.legacy.Move import RigidBodyMC
from sire.search import approx_equal


def test_temperature():
    temp = convert(100, fahrenheit, celsius)

    print("100 F == %f C" % (100*fahrenheit).to(celsius))

    assert approx_equal(temp, (100-32)/1.8)
    assert approx_equal(temp, (100*fahrenheit).to(celsius))

    mc = RigidBodyMC()

    mc.set_temperature( 100 * fahrenheit )

    print(mc.temperature().to(celsius))

    assert approx_equal(mc.temperature().to(fahrenheit), 100)

    k = (4 * pi * 8.854187817e-12 * farad / meter)

    print(k)

    k = 1 / k

    print(k)

    assert approx_equal(k.value(), 332.063710)
