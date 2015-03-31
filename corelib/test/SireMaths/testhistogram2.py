
from Sire.Maths import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot

from nose.tools import assert_almost_equal

rand = RanGenerator()

rand_mean = 20
rand_variance = 5

h1 = Histogram(0.1)
h2 = Histogram(0.1)

nvals = 100000
weight = 1

print("Generating random numbers...")
for i in range(0,nvals):
    h1.accumulate( rand.randNorm(rand_mean, rand_variance), weight )
    h2.accumulate( rand.randNorm(rand_mean+20, rand_variance), weight )

h3 = h1 + h2

def plot(axes, hist, type):
    x = []
    y = []

    for value in list(hist.values()):
        x.append( value.middle() )
        y.append( value.value() )

    axes.plot(x, y, type)

figure = pyplot.figure()
axes = figure.add_subplot(111)

plot(axes, h1, "r-")
plot(axes, h2, "b-")
plot(axes, h3, "g-")

pyplot.show()


