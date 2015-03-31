
from Sire.Soiree import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Units import *

lam = Symbol("lambda")

w = WHAM(lam)

rand = RanGenerator()

coords = []

for i in range(0,100000):
    coords.append( rand.rand() )

umbrella = []

for c in coords:
    umbrella.append( 10 * c * kcal_per_mol )

w.add( coords, umbrella )

umbrella = []

for c in coords:
    umbrella.append( 10 * (1-c) * kcal_per_mol )

print(w.solve( HistogramRange(0,1,0.1) ))

w.add( coords, umbrella )

print(w.solve( HistogramRange(0,1,0.1) ))
