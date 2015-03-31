
from Sire.Maths import *

array = []
for i in range(1,101):
    array.append( i * 0.03 )
    array.append( -i * 0.03 )

array = MultiFixed.fromArray(array)

print(array)

