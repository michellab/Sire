
from Sire.MM import *
from Sire.IO import *
from Sire.DB import *
from Sire.Qt import *

m = ParameterDB()

t = QTime()

t.start()

ProtoMS().read("/home/chris/ProtoMC/trunk/parameter/solvents.ff", m)
ProtoMS().read("/home/chris/ProtoMC/trunk/parameter/opls96.ff", m)
ProtoMS().read("/home/chris/ProtoMC/trunk/parameter/opls96-residues.ff", m)

ms = t.elapsed()
print("Reading parameter files took %d ms" % ms)

t.start()

m.dump("test.db")

ms = t.elapsed()

print("Dumping database took %d ms" % ms)

print("nrelationships = %d, nparameters = %d" % (m.nRelationships(),m.nParameters()))

m2 = ParameterDB()

t.start()

m2.load("test.db")

ms = t.elapsed()

print("Loading database took %d ms" % ms)

t.start()

m2.dump("test2.db")

ms = t.elapsed()

print("Dumping second database took %d ms" % ms)

print("nrelationships = %d, nparameters = %d" % (m2.nRelationships(),m2.nParameters()))

print("test.db and test2.db should be identical - check using diff")

################
## Timings
##
##  icc on cubert:
##  reading parameter files took 1134 ms
##  dumping database took 107 ms
##  loading database took 250 ms
##  dumping database again took 108 ms
##
##  gcc3.4 on cubert
##  reading parameter file took 1117 ms
##  dumping database took 105 ms
##  loading database took 250 ms
##  dumping database again took 108 ms
##
##  gcc4.1 on scruffy
##  reading parameter file took 580 ms
##  dumping database took 42 ms
##  loading database took 116 ms
##  dumping database again took 43 ms
##
##
##  Version 801
##
##  gcc 4.1 on cubert:
##  reading parameter files took 1005 ms
##  dumping database took 154 ms
##  loading database took 266 ms
##  dumping database again took 149 ms
##
##  Version 814 (new ParameterDB version)
##
##  gcc 4.1 on cubert:
##  reading parameter files took 1009 ms
##  dumping database took 157 ms
##  loading database took 273 ms
##  dumping second database took 155 ms

