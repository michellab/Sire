/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008-2013  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "histogram.h"

#include "SireMaths/maths.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMaths;
using namespace SireID;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of HistogramBin
///////////

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const HistogramBin &bin)
{
    ds << bin.minimum() << bin.maximum();
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, HistogramBin &bin)
{
    ds >> bin.minval >> bin.maxval;
    return ds;
}

/** Null constructor */
HistogramBin::HistogramBin() : minval(0), maxval(0)
{}

/** Construct a bin that contains the values that match
    minval <= value < maxval */
HistogramBin::HistogramBin(double min, double max)
             : minval(min), maxval(max)
{
    if (minval > maxval)
    {
        qSwap(minval, maxval);
    }
}

/** Copy constructor */
HistogramBin::HistogramBin(const HistogramBin &other)
             : minval(other.minval), maxval(other.maxval)
{}

/** Copy assignment operator */
HistogramBin& HistogramBin::operator=(const HistogramBin &other)
{
    minval = other.minval;
    maxval = other.maxval;
    return *this;
}

/** Comparison operator */
bool HistogramBin::operator==(const HistogramBin &other) const
{
    return minval == other.minval and maxval == other.maxval;
}

/** Comparison operator */
bool HistogramBin::operator!=(const HistogramBin &other) const
{
    return minval != other.minval or maxval != other.maxval;
}

/** Return the minimum value of the bin */
double HistogramBin::minimum() const
{
    return minval;
}

/** Return the value at the middle of the bin */
double HistogramBin::middle() const
{
    return 0.5*maxval + 0.5*minval;
}

/** Return the maximum value of the bin */
double HistogramBin::maximum() const
{
    return maxval;
}

/** Return a string representation */
QString HistogramBin::toString() const
{
    return QObject::tr("Bin[ %1 <= x < %2 ]")
                    .arg(minval).arg(maxval);
}

///////////
/////////// Implementation of HistogramValue
///////////

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const HistogramValue &value)
{
    ds << value.val << static_cast<const HistogramBin&>(value);
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds,
                                         HistogramValue &value)
{
    ds >> value.val >> static_cast<HistogramBin&>(value);
    return ds;
}

/** Null constructor */
HistogramValue::HistogramValue() : HistogramBin(), val(0)
{}

/** Construct the value for the bin 'bin' equal to 'value' */
HistogramValue::HistogramValue(const HistogramBin &bin, double value)
               : HistogramBin(bin), val(value)
{}

/** Copy constructor */
HistogramValue::HistogramValue(const HistogramValue &other)
               : HistogramBin(other), val(other.val)
{}

/** Destructor */
HistogramValue::~HistogramValue()
{}

/** Copy assignment operator */
HistogramValue& HistogramValue::operator=(const HistogramValue &other)
{
    val = other.val;
    HistogramBin::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool HistogramValue::operator==(const HistogramValue &other) const
{
    return val == other.val and HistogramBin::operator==(other);
}

/** Comparison operator */
bool HistogramValue::operator!=(const HistogramValue &other) const
{
    return val != other.val or HistogramBin::operator!=(other);
}

/** Return the value of the bin */
double HistogramValue::value() const
{
    return val;
}

/** Return a string representation */
QString HistogramValue::toString() const
{
    return QObject::tr("Bin[ %1 <= x < %2 ] == %3")
                    .arg( minimum() ).arg( maximum() ).arg(val);
}

///////////
/////////// Implementation of Histogram
///////////

static const RegisterMetaType<Histogram> r_histogram;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Histogram &histogram)
{
    writeHeader(ds, r_histogram, 2);
    
    SharedDataStream sds(ds);
    
    sds << histogram.binvals << histogram.binwidth
        << histogram.avgval << histogram.avgval2
        << histogram.sum_of_bins
        << static_cast<const Property&>(histogram);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Histogram &histogram)
{
    VersionID v = readHeader(ds, r_histogram);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);

        sds >> histogram.binvals >> histogram.binwidth
            >> histogram.avgval >> histogram.avgval2
            >> histogram.sum_of_bins
            >> static_cast<Property&>(histogram);
    }
    else
        throw version_error( v, "2", r_histogram, CODELOC );
        
    return ds;
}

/** Construct an empty histogram with a bin spacing of 1.0 */
Histogram::Histogram()
          : ConcreteProperty<Histogram,Property>(),
            binwidth(1), avgval(0), avgval2(0), sum_of_bins(0)
{}

/** Construct an empty histogram with specified bin width. Note that
    if the binwidth is less than or equal to zero, then a histogram
    will not be collected, and only the mean and standard deviation
    will be recorded */
Histogram::Histogram(double width)
          : ConcreteProperty<Histogram,Property>(),
            binwidth(width), avgval(0), avgval2(0), sum_of_bins(0)
{
    if (width <= 0)
        binwidth = 0;
    else if (width > 1e20)
        binwidth = 1e20;
}

/** Construct a histogram of specified bin width, and populating it with 
    the passed values (which are all assumed to have weight "1")
    Note that if the binwidth is less than or equal to zero, then a histogram
    will not be collected, and only the mean and standard deviation
    will be recorded */
Histogram::Histogram(double width, const QVector<double> &values)
          : ConcreteProperty<Histogram,Property>(),
            binwidth(width), avgval(0), avgval2(0), sum_of_bins(0)
{
    if (width <= 0)
        binwidth = 0;
    else if (width > 1e20)
        binwidth = 1e20;

    this->operator+=(values);
}

/** Copy constructor */
Histogram::Histogram(const Histogram &other)
          : ConcreteProperty<Histogram,Property>(other),
            binvals(other.binvals), binwidth(other.binwidth),
            avgval(other.avgval), avgval2(other.avgval2),
            sum_of_bins(other.sum_of_bins)
{}

/** Destructor */
Histogram::~Histogram()
{}

/** Copy assignment operator */
Histogram& Histogram::operator=(const Histogram &other)
{
    if (this != &other)
    {
        binvals = other.binvals;
        binwidth = other.binwidth;
        avgval = other.avgval;
        avgval2 = other.avgval2;
        sum_of_bins = other.sum_of_bins;
    }
    
    return *this;
}

const char* Histogram::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Histogram>() );
}

const char* Histogram::what() const
{
    return Histogram::typeName();
}

Histogram* Histogram::clone() const
{
    return new Histogram(*this);
}

/** Comparison operator */
bool Histogram::operator==(const Histogram &other) const
{
    return this == &other or
           (binvals == other.binvals and binwidth == other.binwidth and
            avgval == other.avgval and sum_of_bins == other.sum_of_bins and
            avgval2 == other.avgval2);
}

static qint64 getBin(double x, double binwidth)
{
    if (x == 0)
        return 0;
    else if (x > 0)
        return qint64( x / binwidth );
    else
        return qint64( x / binwidth ) - 1;
}

/** Comparison operator */
bool Histogram::operator!=(const Histogram &other) const
{
    return not Histogram::operator==(other);
}

/** Add the contects of the passed histogram into this histogram */
Histogram& Histogram::operator+=(const Histogram &other)
{
    this->accumulate(other);
    return *this;
}

/** Add the passed value into this histogram */
Histogram& Histogram::operator+=(double value)
{
    this->accumulate(value);
    return *this;
}

/** Add the passed array of values onto this histogram */
Histogram& Histogram::operator+=(const QVector<double> &values)
{
    this->accumulate(values);
    return *this;
}

/** Return the accumulation of the two passed histograms. Note that the
    returned histogram will have a bin width that is equal to the smallest
    bin width of this or other */
Histogram Histogram::operator+(const Histogram &other) const
{
    if (binwidth <= other.binwidth)
    {
        Histogram ret(*this);
        ret += other;
        return ret;
    }
    else
    {
        Histogram ret(other);
        ret += *this;
        return ret;
    }
}

/** Return the histogram that is a copy of this, but on which 'value' has 
    been added */
Histogram Histogram::operator+(double value) const
{
    Histogram ret(*this);
    ret += value;
    return ret;
}

/** Return the histogram that is a copy of this, but on which the passed
    values have been added */
Histogram Histogram::operator+(const QVector<double> &values) const
{
    Histogram ret(*this);
    ret += values;
    return ret;
}

/** Return a string representation of this histogram */
QString Histogram::toString() const
{
    return QObject::tr("Histogram( binWidth() => %1, mean() => %2, sumOfBins() => %3 )")
                .arg(binWidth()).arg(mean()).arg(sumOfBins());
}

/** Return the ith bin in the histogram */
HistogramValue Histogram::operator[](int i) const
{
    i = Index(i).map(binvals.count());
    return values().at(i);
}

/** Return the ith bin in the histogram */
HistogramValue Histogram::at(int i) const
{
    return this->operator[](i);
}

/** Return the number of bins in the histogram */
int Histogram::count() const
{
    return binvals.count();
}

/** Return the number of bins in the histogram */
int Histogram::size() const
{
    return this->count();
}

/** Return the set of all bins and values in the histogram. The bins
    will be returned in numerical order */
QVector<HistogramValue> Histogram::values() const
{
    if (binvals.isEmpty())
        return QVector<HistogramValue>();

    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    QVector<HistogramValue> vals(binvals.count());
    HistogramValue *val = vals.data();
    
    foreach (qint64 bin, bins)
    {
        *val = HistogramValue( HistogramBin(bin*binwidth, (bin+1)*binwidth),
                               binvals[bin] );
        
        val += 1;
    }
    
    return vals;
}

/** Return the idealised normal distribution for the values in the histogram,
    based on the current mean and standard deviation, and the sum of weights */
QVector<HistogramValue> Histogram::normalDistribution() const
{
    if (binvals.isEmpty())
        return QVector<HistogramValue>();
    
    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    QVector<HistogramValue> vals(binvals.count());
    HistogramValue *val = vals.data();
    
    const double avg = this->mean();
    const double stdev = this->standardDeviation();
    const double denom = 1.0 / (2*stdev*stdev);
    double norm = 0;
    int nnorm = 0;

    foreach (qint64 bin, bins)
    {
        double x = (bin+0.5)*binwidth;
        
        if (std::abs( avg - x ) < stdev)
        {
            norm += binvals.value(bin) / std::exp( -pow_2(x-avg) * denom );
            nnorm += 1;
        }
    }

    if (nnorm > 0)
    {
        norm /= nnorm;
    }
    else
    {
        foreach (qint64 bin, bins)
        {
            double x = (bin+0.5)*binwidth;
            norm += binvals.value(bin) / std::exp( -pow_2(x-avg) * denom );
        }
        
        norm /= binvals.count();
    }
    
    foreach (qint64 bin, bins)
    {
        double x = (bin+0.5)*binwidth;
        
        *val = HistogramValue( HistogramBin(bin*binwidth, (bin+1)*binwidth),
                               norm * std::exp( -pow_2(x-avg) * denom ) );
        
        val += 1;
    }

    return vals;
}

/** Accumulate 'value' onto the histogram */
void Histogram::accumulate(double value)
{
    this->accumulate(value, 1.0);
}

/** Accumulate 'value' with the passed 'weight' onto the histogram */
void Histogram::accumulate(double value, double weight)
{
    //we cannot add negative weight to the histogram
    if (weight <= 0)
        return;
    
    //first, calculate the average of the
    if (sum_of_bins == 0)
    {
        avgval = value;
        avgval2 = value*value;
        sum_of_bins = weight;
    }
    else
    {
        const double bigratio = sum_of_bins / (sum_of_bins + weight);
        const double smallratio = 1.0 - bigratio;
    
        avgval = bigratio*avgval + smallratio*value;
        avgval2 = bigratio*avgval2 + smallratio*value*value;
        sum_of_bins += weight;
    }
    
    if (binwidth > 0)
    {
        //now histogram the data
        qint64 bin = getBin(value, binwidth);
    
        binvals.insert(bin, binvals.value(bin,0) + weight);
    }
}

/** Accumulate the passed values onto this histogram */
void Histogram::accumulate(const QVector<double> &values)
{
    foreach (double value, values)
    {
        this->accumulate(value, 1);
    }
}

/** Accumulate the data from the passed histogram onto this histogram */
void Histogram::accumulate(const Histogram &other)
{
    if (this->binWidth() > 0 and other.binWidth() > 0)
    {
        Histogram resized = other.resize( this->binWidth() );
    
        for (QHash<qint64,double>::const_iterator it = resized.binvals.constBegin();
             it != resized.binvals.constEnd();
             ++it)
        {
            this->accumulate( (it.key() + 0.5)*resized.binwidth, it.value() );
        }
    }
    else
    {
        //one of the histograms has a binwidth of 0. This destroys the
        //histogram data, leaving only the mean and standard deviation
        binvals.clear();
        binwidth = 0;

        const double bigratio = sum_of_bins / (sum_of_bins + other.sum_of_bins);
        const double smallratio = 1.0 - bigratio;
    
        avgval = bigratio*avgval + smallratio*other.avgval;
        avgval2 = bigratio*avgval2 + smallratio*other.avgval2;
        sum_of_bins += other.sum_of_bins;
    }
}

/** Add 'value' onto the histogram */
void Histogram::add(double value)
{
    this->accumulate(value);
}

/** Add the passed values on this histogram */
void Histogram::add(const QVector<double> &values)
{
    this->accumulate(values);
}

/** Add 'value' with the passed 'weight' onto the histogram */
void Histogram::add(double value, double weight)
{
    this->accumulate(value, weight);
}

/** Add the passed histogram onto this histogram. This will match the
    bin width of the passed histogram to this histogram */
void Histogram::add(const Histogram &other)
{
    this->accumulate(other);
}

/** Return the sum of the weights over all of the bins */
double Histogram::sumOfBins() const
{
    return sum_of_bins;
}

/** Return the width of the bins */
double Histogram::binWidth() const
{
    return binwidth;
}

/** Return the mean average of all values added to the histogram. This
    is calculated exactly from the added data */
double Histogram::mean() const
{
    return avgval;
}

/** Return the mean of the square values */
double Histogram::meanOfSquares() const
{
    return avgval2;
}

/** Return the standard deviation of all values added to the histogram.
    This is calculated exactly from the added data */
double Histogram::standardDeviation() const
{
    if (binvals.isEmpty())
        return 0;

    return std::sqrt(avgval2 - (avgval*avgval));
}

/** Return the standard error of the mean (standard deviation 
    divided by the square root of the number of samples) */
double Histogram::standardError() const
{
    if (sum_of_bins == 0)
        return 0;
    else
        return standardDeviation() / std::sqrt(sum_of_bins);
}

/** Return the students t-value for the passed confidence level
    for the passed number of samples for the mean */
double Histogram::tValue(int nsamples, double level)
{
    // data copied from wikipedia - "http://en.wikipedia.org/wiki/Student's_t-distribution"
    // see t_values and make_tvalues.py in the same directory as this source
    int nlevels = 11;
    
    double levels[11] = { 50.0, 60.0, 70.0, 80.0, 90.0, 95.0, 98.0, 99.0, 99.5, 99.8, 99.9 };
    
    int ncounts = 37;
    
    int counts[37] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40, 50, 60, 80, 100, 120, 1000 };
    
    double values[37][11] = {
        { 1.0, 1.376, 1.963, 3.078, 6.314, 12.71, 31.82, 63.66, 127.3, 318.3, 636.6 },
        { 0.816, 1.061, 1.386, 1.886, 2.92, 4.303, 6.965, 9.925, 14.09, 22.33, 31.6 },
        { 0.765, 0.978, 1.25, 1.638, 2.353, 3.182, 4.541, 5.841, 7.453, 10.21, 12.92 },
        { 0.741, 0.941, 1.19, 1.533, 2.132, 2.776, 3.747, 4.604, 5.598, 7.173, 8.61 },
        { 0.727, 0.92, 1.156, 1.476, 2.015, 2.571, 3.365, 4.032, 4.773, 5.893, 6.869 },
        { 0.718, 0.906, 1.134, 1.44, 1.943, 2.447, 3.143, 3.707, 4.317, 5.208, 5.959 },
        { 0.711, 0.896, 1.119, 1.415, 1.895, 2.365, 2.998, 3.499, 4.029, 4.785, 5.408 },
        { 0.706, 0.889, 1.108, 1.397, 1.86, 2.306, 2.896, 3.355, 3.833, 4.501, 5.041 },
        { 0.703, 0.883, 1.1, 1.383, 1.833, 2.262, 2.821, 3.25, 3.69, 4.297, 4.781 },
        { 0.7, 0.879, 1.093, 1.372, 1.812, 2.228, 2.764, 3.169, 3.581, 4.144, 4.587 },
        { 0.697, 0.876, 1.088, 1.363, 1.796, 2.201, 2.718, 3.106, 3.497, 4.025, 4.437 },
        { 0.695, 0.873, 1.083, 1.356, 1.782, 2.179, 2.681, 3.055, 3.428, 3.93, 4.318 },
        { 0.694, 0.87, 1.079, 1.35, 1.771, 2.16, 2.65, 3.012, 3.372, 3.852, 4.221 },
        { 0.692, 0.868, 1.076, 1.345, 1.761, 2.145, 2.624, 2.977, 3.326, 3.787, 4.14 },
        { 0.691, 0.866, 1.074, 1.341, 1.753, 2.131, 2.602, 2.947, 3.286, 3.733, 4.073 },
        { 0.69, 0.865, 1.071, 1.337, 1.746, 2.12, 2.583, 2.921, 3.252, 3.686, 4.015 },
        { 0.689, 0.863, 1.069, 1.333, 1.74, 2.11, 2.567, 2.898, 3.222, 3.646, 3.965 },
        { 0.688, 0.862, 1.067, 1.33, 1.734, 2.101, 2.552, 2.878, 3.197, 3.61, 3.922 },
        { 0.688, 0.861, 1.066, 1.328, 1.729, 2.093, 2.539, 2.861, 3.174, 3.579, 3.883 },
        { 0.687, 0.86, 1.064, 1.325, 1.725, 2.086, 2.528, 2.845, 3.153, 3.552, 3.85 },
        { 0.686, 0.859, 1.063, 1.323, 1.721, 2.08, 2.518, 2.831, 3.135, 3.527, 3.819 },
        { 0.686, 0.858, 1.061, 1.321, 1.717, 2.074, 2.508, 2.819, 3.119, 3.505, 3.792 },
        { 0.685, 0.858, 1.06, 1.319, 1.714, 2.069, 2.5, 2.807, 3.104, 3.485, 3.767 },
        { 0.685, 0.857, 1.059, 1.318, 1.711, 2.064, 2.492, 2.797, 3.091, 3.467, 3.745 },
        { 0.684, 0.856, 1.058, 1.316, 1.708, 2.06, 2.485, 2.787, 3.078, 3.45, 3.725 },
        { 0.684, 0.856, 1.058, 1.315, 1.706, 2.056, 2.479, 2.779, 3.067, 3.435, 3.707 },
        { 0.684, 0.855, 1.057, 1.314, 1.703, 2.052, 2.473, 2.771, 3.057, 3.421, 3.69 },
        { 0.683, 0.855, 1.056, 1.313, 1.701, 2.048, 2.467, 2.763, 3.047, 3.408, 3.674 },
        { 0.683, 0.854, 1.055, 1.311, 1.699, 2.045, 2.462, 2.756, 3.038, 3.396, 3.659 },
        { 0.683, 0.854, 1.055, 1.31, 1.697, 2.042, 2.457, 2.75, 3.03, 3.385, 3.646 },
        { 0.681, 0.851, 1.05, 1.303, 1.684, 2.021, 2.423, 2.704, 2.971, 3.307, 3.551 },
        { 0.679, 0.849, 1.047, 1.299, 1.676, 2.009, 2.403, 2.678, 2.937, 3.261, 3.496 },
        { 0.679, 0.848, 1.045, 1.296, 1.671, 2.0, 2.39, 2.66, 2.915, 3.232, 3.46 },
        { 0.678, 0.846, 1.043, 1.292, 1.664, 1.99, 2.374, 2.639, 2.887, 3.195, 3.416 },
        { 0.677, 0.845, 1.042, 1.29, 1.66, 1.984, 2.364, 2.626, 2.871, 3.174, 3.39 },
        { 0.677, 0.845, 1.041, 1.289, 1.658, 1.98, 2.358, 2.617, 2.86, 3.16, 3.373 },
        { 0.674, 0.842, 1.036, 1.282, 1.645, 1.96, 2.326, 2.576, 2.807, 3.09, 3.291 }
    };
    
    //first find the index for the number of samples
    int sample_idx = ncounts - 1;
    
    for (int i=0; i<ncounts; ++i)
    {
        if (nsamples < counts[i])
        {
            sample_idx = i-1;
            break;
        }
    }
    
    if (sample_idx < 0)
        sample_idx = 0;
    
    //now find the confidence level
    int level_idx = nlevels - 1;
    
    for (int i=0; i<nlevels; ++i)
    {
        if (level < levels[i])
        {
            level_idx = i - 1;
            break;
        }
    }
    
    if (level_idx < 0)
        level_idx = 0;
    
    return values[sample_idx][level_idx];
}

/** Return the students t-value for the passed confidence level
    for the number of samples in the histogram */
double Histogram::tValue(double level) const
{
    return tValue( int(sum_of_bins+0.5), level );
}

/** Return the standard error calculated to the passed level 
    (66, 90, 95 or 99%) */
double Histogram::standardError(double level) const
{
    return tValue(level) * standardError();
}

/** Return the skew of the data. This is estimated based on the histogram
    of the data */
double Histogram::skew() const
{
    double avg = mean();
    double stdev = standardDeviation();
    double denom = 1.0 / (sum_of_bins*pow_3(stdev));
    
    double skew = 0;
    
    for (QHash<qint64,double>::const_iterator it = binvals.constBegin();
         it != binvals.constEnd();
         ++it)
    {
        skew += denom * it.value() * pow_3((it.key()+0.5) * binwidth - avg);
    }
    
    return skew;
}

/** Return the excess kirtosis of the data. This is estimated based on the histogram
    of the data (this is the kirtosis minus 3, so that the normal distribution
    has a kirtosis of 0) */
double Histogram::kirtosis() const
{
    double avg = mean();
    double stdev = standardDeviation();
    double denom = 1.0 / (sum_of_bins*pow_4(stdev));
    
    double kirt = 0;
    
    for (QHash<qint64,double>::const_iterator it = binvals.constBegin();
         it != binvals.constEnd();
         ++it)
    {
        kirt += denom * it.value() * pow_4( (it.key() + 0.5)*binwidth - avg );
    }
    
    return kirt - 3;
}

/** Return the median of all values added to the histogram. This is
    estimated based on the actual histogram of added data */
double Histogram::median() const
{
    if (binvals.isEmpty())
        return 0;

    double sum = 0;
    const double half_full = 0.5 * sumOfBins();
    
    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    foreach (qint64 bin, bins)
    {
        sum += binvals[bin];
        
        if (sum > half_full)
        {
            //by how much have we gone over...
            double amount = (sum - half_full) / binvals[bin];
            
            return (bin+amount)*binwidth;
        }
    }
    
    throw SireError::program_bug( QObject::tr(
            "It should not be possible to reach here...!"), CODELOC );
    
    return 0;
}

/** Return the mode of all values added to the histogram. This is 
    estimated based on the actual histogram of added data */
double Histogram::mode() const
{
    if (binvals.isEmpty())
        return 0;
    
    double maxval = 0;
    qint64 maxbin = 0;
    
    for (QHash<qint64,double>::const_iterator it = binvals.constBegin();
         it != binvals.constEnd();
         ++it)
    {
        if (it.value() > maxval)
        {
            maxval = it.value();
            maxbin = it.key();
        }
    }
    
    return (maxbin+0.5) * binwidth;
}

/** Return the highest value in the histogram */
double Histogram::maximumValue() const
{
    if (binvals.isEmpty())
        return 0;
    
    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    return (bins.last() + 1) * binwidth;
}

/** Return the lowest values in the histogram */
double Histogram::minimumValue() const
{
    if (binvals.isEmpty())
        return 0;
    
    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    return bins.first() * binwidth;
}

/** Return the range for the data in the histogram */
double Histogram::range() const
{
    if (binvals.isEmpty())
        return 0;
    
    QList<qint64> bins = binvals.keys();
    qSort(bins);
    
    return (bins.last() - bins.first() + 1) * binwidth;
}

/** Return a normalised version of this histogram. The histogram
    is normalised so that the sum under the curve is 1 (e.g.
    sum_of_bins * bin_width is 1) */
Histogram Histogram::normalise() const
{
    if (binvals.isEmpty())
        return Histogram();

    else if (binwidth / sum_of_bins == 1)
        return *this;

    Histogram ret(*this);
    
    foreach (qint64 bin, ret.binvals.keys())
    {
        ret.binvals.insert(bin, ret.binvals.value(bin) / (binwidth*sum_of_bins));
    }
    
    ret.sum_of_bins = binwidth;
    
    return ret;
}

/** Return a resized copy of this histogram with the passed new binwidth */
Histogram Histogram::resize(double width) const
{
    if (width == binwidth)
        return *this;

    Histogram ret;
    ret.binwidth = width;
    ret.avgval = avgval;
    ret.avgval2 = avgval2;
    ret.sum_of_bins = sum_of_bins;

    if (width <= 0)
    {
        ret.binwidth = 0;
        return ret;
    }

    if (width > 1e20)
    {
        ret.binwidth = 1e20;
    }
    
    for (QHash<qint64,double>::const_iterator it = binvals.constBegin();
         it != binvals.constEnd();
         ++it)
    {
        double weight = it.value();
    
        double old_minval = it.key() * binwidth;
        double old_maxval = old_minval + binwidth;

        qint64 bin = getBin(old_minval, width);

        while (weight > 0)
        {
            double new_maxval = (bin+1)*width;
        
            if (new_maxval < old_maxval)
            {
                double partial_weight = weight * (new_maxval - old_minval) / binwidth;
                ret.binvals.insert( bin, ret.binvals.value(bin,0) + partial_weight );
                weight -= partial_weight;
                old_minval = new_maxval;
            }
            else
            {
                ret.binvals.insert( bin, ret.binvals.value(bin,0) + weight );
                weight = 0;
            }
            
            bin += 1;
        }
    }
    
    if (width < binwidth)
    {
        //the new histogram has a higher resolution, so will need to be smoothed
        QHash<qint64,double> smoothed = ret.binvals;
        QList<qint64> bins = smoothed.keys();
        qSort(bins);
        
        for (int i=1; i<bins.count()-1; ++i)
        {
            smoothed[ bins[i] ] = 0.25*ret.binvals[bins[i-1]] + 0.5 * ret.binvals[bins[i]]
                                      + 0.25*ret.binvals[bins[i+1]];
        }
        
        smoothed[bins[0]] = 0.75*ret.binvals[bins[0]] + 0.25*ret.binvals[bins[1]];
        
        if (bins.count() > 1)
        {
            smoothed[bins[bins.count()-1]] = 0.75*ret.binvals[bins[bins.count()-1]] +
                                0.25 * ret.binvals[bins[bins.count()-2]];
        }
        
        ret.binvals = smoothed;
    }
    
    return ret;
}
