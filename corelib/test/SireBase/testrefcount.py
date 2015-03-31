
from Sire.Base import *

mieow = StringProperty("mieow")
p = Properties()
p.setProperty("cat", mieow)

print(p)

print("Deleting 'mieow' - should not crash!")
mieow = 0

print("Printing 'p' - should not crash!")
print(p)

print("Deleting 'p' - should not crash!")
p = 0

print("Everything worked")
