
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot

from Sire.Maths import Histogram, HistogramValue
from Sire.Analysis import DataPoint

def _plotHistogram(histogram):
    x = []
    y = []
    for value in list(histogram.values()):
        x.append(value.minimum())
        y.append(0)
        x.append(value.minimum())
        y.append(value.value())
        x.append(value.maximum())
        y.append(value.value())
        x.append(value.maximum())
        y.append(0)

    pyplot.errorbar(x, y)

def _plotDataPoint(point):
    x = [point.x()]
    y = [point.y()]

    if point.hasErrorRange():
        xerr = [point.xMinError()]
        yerr = [point.yMinError()]
        pyplot.errorbar(x, y, yerr, xerr)
        xerr = [point.xMaxError()]
        yerr = [point.yMaxError()]
        pyplot.errorbar(x, y, yerr, xerr)
    elif point.hasError():
        xerr = [point.xError()]
        yerr = [point.yError()]
        pyplot.errorbar(x, y, yerr, xerr)
    else:
        pyplot.errorbar(x, y)

def _tryPlotHistValues(points):
    x = []
    y = []

    for point in points:
        if not isinstance(point, HistogramValue):
            return False

        x.append(point.middle())
        y.append(point.value())

    pyplot.errorbar(x,y)

    return True

def _tryPlotDataPoints(points):
    x = []
    y = []
    xminerr = []
    xmaxerr = []
    yminerr = []
    ymaxerr = []

    has_error_range = False
    has_error = False

    for point in points:
        if not isinstance(point, DataPoint):
            return False

        x.append(point.x())
        y.append(point.y())

        if point.hasErrorRange():
            has_error_range = True
            xminerr.append( point.xMinError() )
            xmaxerr.append( point.xMaxError() )
            yminerr.append( point.yMinError() )
            ymaxerr.append( point.yMaxError() )
        elif point.hasError():
            has_error = True
            xminerr.append( point.xError() )
            xmaxerr.append( point.xError() )
            yminerr.append( point.yError() )
            ymaxerr.append( point.yError() )
        else:
            xminerr.append(0)
            xmaxerr.append(0)
            yminerr.append(0)
            ymaxerr.append(0)

    if has_error_range:
        pyplot.errorbar(x, y, ymaxerr, xmaxerr)
        pyplot.errorbar(x, y, yminerr, xminerr)
    elif has_error:
        pyplot.errorbar(x, y, ymaxerr, xmaxerr)
    else:
        pyplot.errorbar(x, y)
    
    return True

def _plot(graph):
    if isinstance(graph, Histogram):
        _plotHistogram(graph)
    
    elif isinstance(graph, DataPoint):
        _plotDataPoint(graph)
    
    elif not _tryPlotDataPoints(graph):
        if not _tryPlotHistValues(graph):
            for item in graph:
                _plot(item)

def plot(graph):
    _plot(graph)
    pyplot.show()

