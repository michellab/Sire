
from Sire.Maths import *

import matplotlib
matplotlib.use('TkAgg')

from matplotlib import pyplot

from nose.tools import assert_almost_equal

rand = RanGenerator()

rand_mean = 20
rand_variance = 5

h = Histogram(0.2)

assert h.binWidth() == 0.2

nvals = 100000
weight = 0.6

print("Generating random numbers...")
for i in range(0,nvals):
    h.accumulate( rand.randNorm(rand_mean, rand_variance), weight )

assert_almost_equal( h.sumOfBins(), nvals * weight, 1 )
assert_almost_equal( h.mean(), rand_mean, 1 )
assert_almost_equal( h.standardDeviation(), rand_variance, 1 )

print(h)
print(("MEAN", h.mean(), "vs.", rand_mean))
print(("MEDIAN", h.median()))
print(("MODE", h.mode()))
print(("RANGE", h.range()))
print(("MINVAL", h.minimumValue()))
print(("MAXVAL", h.maximumValue()))
print(("STDEV", h.standardDeviation(), "vs.", rand_variance))
print(("SKEW", h.skew(), "vs.", 0.0))
print(("KIRTOSIS", h.kirtosis(), "vs.", 0.0))

vals = list(h.normalise().values())
gaus = h.normalise().normalDistribution()

assert len(vals) == len(gaus)

x = []
y1 = []
y2 = []

for i in range(0,len(vals)):
    x.append( vals[i].middle() )
    y1.append( vals[i].value() )
    y2.append( gaus[i].value() )

figure = pyplot.figure()
axes = figure.add_subplot(111)
axes.plot(x, y1, "b-")
axes.plot(x, y2, "r-")

x = []
y = []
vals = list(h.resize(0.4).normalise().values())

for i in range(0,len(vals)):
    x.append( vals[i].middle() )
    y.append( vals[i].value() )

axes.plot(x, y, "g-")

x = []
y = []
vals = list(h.resize(0.1).normalise().values())

for i in range(0,len(vals)):
    x.append( vals[i].middle() )
    y.append( vals[i].value() )

axes.plot(x, y, "y-")

pyplot.show()
