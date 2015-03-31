
from Sire.Units import *
from Sire.Move import *

from Sire.Maths import pi

def close( val0, val1 ):
   
   if not (abs(float(val1)-float(val0)) < 0.000001):
        print("FAIL! %f != %f" % (val0, val1))
        return False
   else:
        return True

temp = convert(100, fahrenheit, celsius)

print("100 F == %f C" % (100*fahrenheit).to(celsius))

assert( close(temp, (100-32)/1.8) )
assert( close(temp, (100*fahrenheit).to(celsius)) )

mc = RigidBodyMC()

mc.setTemperature( 100 * fahrenheit )

print(mc.temperature().to(celsius))

assert( close(mc.temperature().to(fahrenheit), 100) )

k = (4 * pi * 8.854187817e-12 * farad / meter)
print(k)

k = 1 / k
print(k)

assert( close(k.value(), 332.063710) )

