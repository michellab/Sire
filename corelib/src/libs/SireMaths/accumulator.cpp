/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QMutex>

#include <cmath>

#include "accumulator.h"
#include "histogram.h"

#include "SireMaths/maths.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMaths;
using namespace SireBase;
using namespace SireStream;

/////////
///////// Implementation of Accumulator
/////////

static const RegisterMetaType<Accumulator> r_accum( MAGIC_ONLY,
                                                    Accumulator::typeName() );
                                                    
/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Accumulator &accum)
{
    writeHeader(ds, r_accum, 1);
    
    ds << accum.nvalues << static_cast<const Property&>(accum);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Accumulator &accum)
{
    VersionID v = readHeader(ds, r_accum);
    
    if (v == 1)
    {
        ds >> accum.nvalues >> static_cast<Property&>(accum);
    }
    else
        throw version_error(v, "1", r_accum, CODELOC);
        
    return ds;
}

/** Constructor */
Accumulator::Accumulator() : Property(), nvalues(0)
{}

/** Copy constructor */
Accumulator::Accumulator(const Accumulator &other)
            : Property(other), nvalues(other.nvalues)
{}

/** Destructor */
Accumulator::~Accumulator()
{}

/** Return the number of values that have been sampled */
int Accumulator::nSamples() const
{
    return int(nvalues);
}

/** Accumulate the value 'value' onto the sample */
void Accumulator::accumulate(double)
{
    ++nvalues;
}

/** Accumulate many values */
void Accumulator::accumulate(const QVector<double> &values)
{
    foreach (double value, values)
    {
        this->accumulate(value);
    }
}

/** Accumulate many values */
void Accumulator::accumulate(const QList<double> &values)
{
    foreach (double value, values)
    {
        this->accumulate(value);
    }
}

/** Completely clear the statistics in this accumulator */
void Accumulator::clear()
{
    nvalues = 0;
}

/** Internal function used to increment the number of samples */
void Accumulator::add(int nsteps)
{
    nvalues += nsteps;
}

/** Internal copy assignment operator */
Accumulator& Accumulator::operator=(const Accumulator &other)
{
    nvalues = other.nvalues;
    Property::operator=(other);
    return *this;
}

/** Internal comparison operator */
bool Accumulator::operator==(const Accumulator &other) const
{
    return nvalues == other.nvalues;
}

/** Internal comparison operator */
bool Accumulator::operator!=(const Accumulator &other) const
{
    return nvalues != other.nvalues;
}

/////////
///////// Implementation of NullAccumulator
/////////

static const RegisterMetaType<NullAccumulator> r_null;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const NullAccumulator &null)
{
    writeHeader(ds, r_null, 1);
    
    ds << static_cast<const Accumulator&>(null);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, NullAccumulator &null)
{
    VersionID v = readHeader(ds, r_null);
    
    if (v == 1)
    {
        ds >> static_cast<Accumulator&>(null);
    }
    else
        throw version_error(v, "1", r_null, CODELOC);
        
    return ds;
}

/** Construct an empty average */
NullAccumulator::NullAccumulator() 
                : ConcreteProperty<NullAccumulator,Accumulator>()
{}

/** Copy constructor */
NullAccumulator::NullAccumulator(const NullAccumulator &other)
        : ConcreteProperty<NullAccumulator,Accumulator>(other)
{}

/** Destructor */
NullAccumulator::~NullAccumulator()
{}

/** Copy assignment operator */
NullAccumulator& NullAccumulator::operator=(const NullAccumulator &other)
{
    Accumulator::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullAccumulator::operator==(const NullAccumulator &other) const
{
    return true;
}

/** Comparison operator */
bool NullAccumulator::operator!=(const NullAccumulator &other) const
{
    return false;
}

/** Accumulate the passed value onto the average */
void NullAccumulator::accumulate(double)
{}

const char* NullAccumulator::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullAccumulator>() );
}

const NullAccumulator& Accumulator::null()
{
    return *(create_shared_null<NullAccumulator>());
}

/////////
///////// Implementation of Average
/////////

static const RegisterMetaType<Average> r_avg;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Average &avg)
{
    writeHeader(ds, r_avg, 1);
    
    ds << avg.avgval << static_cast<const Accumulator&>(avg);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Average &avg)
{
    VersionID v = readHeader(ds, r_avg);
    
    if (v == 1)
    {
        ds >> avg.avgval >> static_cast<Accumulator&>(avg);
    }
    else
        throw version_error(v, "1", r_avg, CODELOC);
        
    return ds;
}

/** Construct an empty average */
Average::Average() : ConcreteProperty<Average,Accumulator>(), avgval(0)
{}

/** Copy constructor */
Average::Average(const Average &other)
        : ConcreteProperty<Average,Accumulator>(other), avgval(other.avgval)
{}

/** Destructor */
Average::~Average()
{}

/** Copy assignment operator */
Average& Average::operator=(const Average &other)
{
    if (this != &other)
    {
        avgval = other.avgval;
        Accumulator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Average::operator==(const Average &other) const
{
    return avgval == other.avgval and
           Accumulator::operator==(other);
}

/** Comparison operator */
bool Average::operator!=(const Average &other) const
{
    return not this->operator==(other);
}

/** Self-addition operator */
Average& Average::operator+=(const Average &other)
{
    if (nSamples() == 0)
        this->operator=(other);
    
    else if (other.nSamples() > 0)
    {
        double nsteps = nSamples() + other.nSamples();
        
        double my_ratio = nSamples() / nsteps;
        double other_ratio = other.nSamples() / nsteps;
        
        avgval = avgval * my_ratio + other.avgval * other_ratio;
        
        Accumulator::add(other.nSamples());
    }
    
    return *this;
}

/** Addition operator */
Average Average::operator+(const Average &other) const
{
    Average ret(*this);
    ret += other;
    return ret;
}

/** Completely clear the statistics in this accumulator */
void Average::clear()
{
    avgval = 0;
    Accumulator::clear();
}

/** Accumulate the passed value onto the average */
void Average::accumulate(double value)
{
    double nsteps = this->nSamples() + 1;
    
    //calculate the average as
    // average = ((n-1)/n) * average + (1/n) * value
    
    double big_ratio = (nsteps - 1) / nsteps;
    double small_ratio = 1.0 / nsteps;
    
    avgval = (big_ratio * avgval) + (small_ratio * value);
    Accumulator::accumulate(value);
}

/** Return the average value */
double Average::average() const
{
    return avgval;
}

/** Allow automatic casting to a double to retrieve the average value */
Average::operator double() const
{
    return this->average();
}

const char* Average::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Average>() );
}

/////////
///////// Implementation of AverageAndStddev
/////////

static const RegisterMetaType<AverageAndStddev> r_avgstddev;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, 
                                         const AverageAndStddev &avgstddev)
{
    writeHeader(ds, r_avgstddev, 1);
    
    ds << avgstddev.avgval2 << static_cast<const Average&>(avgstddev);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds,
                                         AverageAndStddev &avgstddev)
{
    VersionID v = readHeader(ds, r_avgstddev);
    
    if (v == 1)
    {
        ds >> avgstddev.avgval2 >> static_cast<Average&>(avgstddev);
    }
    else
        throw version_error( v, "1", r_avgstddev, CODELOC );
        
    return ds;
}

/** Construct an empty average */
AverageAndStddev::AverageAndStddev()
                 : ConcreteProperty<AverageAndStddev,Average>(),
                   avgval2(0)
{}

/** Copy constructor */
AverageAndStddev::AverageAndStddev(const AverageAndStddev &other)
                 : ConcreteProperty<AverageAndStddev,Average>(other),
                   avgval2(other.avgval2)
{}

/** Destructor */
AverageAndStddev::~AverageAndStddev()
{}

/** Copy assignment operator */
AverageAndStddev& AverageAndStddev::operator=(const AverageAndStddev &other)
{
    if (this != &other)
    {
        avgval2 = other.avgval2;
        Average::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool AverageAndStddev::operator==(const AverageAndStddev &other) const
{
    return avgval2 == other.avgval2 and
           Average::operator==(other);
}

/** Comparison operator */
bool AverageAndStddev::operator!=(const AverageAndStddev &other) const
{
    return not this->operator==(other);
}

/** Self-addition operator */
AverageAndStddev& AverageAndStddev::operator+=(const AverageAndStddev &other)
{
    if (nSamples() == 0)
        this->operator=(other);

    else if (other.nSamples() > 0)
    {
        Average::operator+=(other);
        
        double nsteps = nSamples() + other.nSamples();
        
        double my_ratio = nSamples() / nsteps;
        double other_ratio = other.nSamples() / nsteps;
        
        avgval2 = avgval2 * my_ratio + other.avgval2 * other_ratio;
    }
    
    return *this;
}

/** Addition operator */
AverageAndStddev AverageAndStddev::operator+(const AverageAndStddev &other) const
{
    AverageAndStddev ret(*this);
    ret += other;
    return ret;
}

/** Completely clear the statistics in this accumulator */
void AverageAndStddev::clear()
{
    avgval2 = 0;
    Average::clear();
}

/** Accumulate the average and standard deviation */
void AverageAndStddev::accumulate(double value)
{
    double nsteps = this->nSamples() + 1;
    
    //calculate the average of the squares as
    // average2 = ((n-1)/n) * average2 + (1/n) * value * value
    
    double big_ratio = (nsteps - 1) / nsteps;
    double small_ratio = 1.0 / nsteps;
    
    avgval2 = (big_ratio * avgval2) + (small_ratio * value * value);

    Average::accumulate(value);
}

/** Return the standard deviation of the average
    (calculated as the sqrt of the mean of the squares minus
     the square of the mean) */
double AverageAndStddev::stddev() const
{
    return std::sqrt( avgval2 - pow_2(this->average()) );
}

/** Return the standard deviation of the average
    (calculated as the sqrt of the mean of the squares minus
     the square of the mean) */
double AverageAndStddev::standardDeviation() const
{
    return this->stddev();
}

/** Return the standard error on the average */
double AverageAndStddev::standardError() const
{
    if (nSamples() == 0)
        return 0;
    else
        return standardDeviation() / sqrt(nSamples());
}

/** Return the standard error calculated to the passed level 
    (66, 90, 95 or 99%) */
double AverageAndStddev::standardError(int level) const
{
    return Histogram::tValue(nSamples(), level) * standardError();
}

/** Return the mean average of the squares */
double AverageAndStddev::meanOfSquares() const
{
    return avgval2;
}

const char* AverageAndStddev::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AverageAndStddev>() );
}

/////////
///////// Implementation of ExpAverage
/////////

static const RegisterMetaType<ExpAverage> r_expavg;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const ExpAverage &expavg)
{
    writeHeader(ds, r_expavg, 2);
    
    ds << expavg.avgval << expavg.avgval2 << expavg.sclfac
       << static_cast<const Accumulator&>(expavg);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, ExpAverage &expavg)
{
    VersionID v = readHeader(ds, r_expavg);
    
    if (v == 2)
    {
        ds >> expavg.avgval >> expavg.avgval2 >> expavg.sclfac
           >> static_cast<Accumulator&>(expavg);
    }
    else if (v == 1)
    {
        ds >> expavg.avgval >> expavg.sclfac
           >> static_cast<Accumulator&>(expavg);
        
        expavg.avgval2 = expavg.avgval * expavg.avgval;;
    }
    else
        throw version_error(v, "1,2", r_expavg, CODELOC);
        
    return ds;
}

/** Construct an empty average using the passed scale factor.

    \throw SireError::invalid_arg
*/
ExpAverage::ExpAverage(double scale_factor) 
           : ConcreteProperty<ExpAverage,Accumulator>(), 
             avgval(0), sclfac(scale_factor)
{}

/** Copy constructor */
ExpAverage::ExpAverage(const ExpAverage &other)
        : ConcreteProperty<ExpAverage,Accumulator>(other), 
          avgval(other.avgval), sclfac(other.sclfac)
{}

/** Destructor */
ExpAverage::~ExpAverage()
{}

/** Copy assignment operator */
ExpAverage& ExpAverage::operator=(const ExpAverage &other)
{
    if (this != &other)
    {
        avgval = other.avgval;
        sclfac = other.sclfac;
        Accumulator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool ExpAverage::operator==(const ExpAverage &other) const
{
    return avgval == other.avgval and
           sclfac == other.sclfac and
           Accumulator::operator==(other);
}

/** Comparison operator */
bool ExpAverage::operator!=(const ExpAverage &other) const
{
    return not this->operator==(other);
}

/** Self-addition operator */
ExpAverage& ExpAverage::operator+=(const ExpAverage &other)
{
    if (sclfac != other.sclfac)
        throw SireError::incompatible_error( QObject::tr(
                "Cannot add together these two ExpAverage objects as their scale "
                "factors are different. %1 vs. %2.")
                    .arg(sclfac).arg(other.sclfac), CODELOC );
    
    if (nSamples() == 0)
        this->operator=(other);
    
    else if (other.nSamples() > 0)
    {
        double nsteps = nSamples() + other.nSamples();
        
        double my_ratio = nSamples() / nsteps;
        double other_ratio = other.nSamples() / nsteps;
        
        avgval = avgval * my_ratio + other.avgval * other_ratio;

        Accumulator::add(other.nSamples());
    }
    
    return *this;
}

/** Addition operator */
ExpAverage ExpAverage::operator+(const ExpAverage &other) const
{
    ExpAverage ret(*this);
    ret += other;
    return ret;
}

/** Completely clear the statistics in this accumulator */
void ExpAverage::clear()
{
    avgval = 0;
    Accumulator::clear();
}

/** Accumulate the passed value onto the average */
void ExpAverage::accumulate(double value)
{
    double expvalue = std::exp( sclfac * value );

    double nsteps = this->nSamples() + 1;
    
    //calculate the average as
    // average = ((n-1)/n) * average + (1/n) * expvalue
    
    double big_ratio = (nsteps - 1) / nsteps;
    double small_ratio = 1.0 / nsteps;
    
    avgval = (big_ratio * avgval) + (small_ratio * expvalue);
    avgval2 = (big_ratio * avgval) + (small_ratio * expvalue*expvalue);
    
    Accumulator::accumulate(value);
}

/** Return the average value */
double ExpAverage::average() const
{
    return ( 1.0 / sclfac ) * std::log(avgval);
}

/** Return the average of the squared value */
double ExpAverage::average2() const
{
    return ( 1.0 / sclfac ) * std::log(std::sqrt(avgval2));
}

/** Allow automatic casting to a double to retrieve the average value */
ExpAverage::operator double() const
{
    return this->average();
}

/** Internal function used to return the scale factor */
double ExpAverage::scaleFactor() const
{
    return sclfac;
}

const char* ExpAverage::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ExpAverage>() );
}

/////////
///////// Implementation of Median
/////////

static const RegisterMetaType<Median> r_median;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Median &median)
{
    writeHeader(ds, r_median, 1);
    
    ds << median.minval << median.maxval
       << static_cast<const Accumulator&>(median);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Median &median)
{
    VersionID v = readHeader(ds, r_median);
    
    if (v == 1)
    {
        ds >> median.minval >> median.maxval
           >> static_cast<Accumulator&>(median);
    }
    else
        throw version_error(v, "1", r_median, CODELOC);
        
    return ds;
}

/** Construct an empty average */
Median::Median() 
       : ConcreteProperty<Median,Accumulator>(), 
         minval(  std::numeric_limits<double>::max() ),
         maxval( -std::numeric_limits<double>::max() )
{}

/** Copy constructor */
Median::Median(const Median &other)
        : ConcreteProperty<Median,Accumulator>(other), 
          minval(other.minval), maxval(other.maxval)
{}

/** Destructor */
Median::~Median()
{}

/** Copy assignment operator */
Median& Median::operator=(const Median &other)
{
    if (this != &other)
    {
        minval = other.minval;
        maxval = other.maxval;
        Accumulator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Median::operator==(const Median &other) const
{
    return minval == other.minval and
           maxval == other.maxval and
           Accumulator::operator==(other);
}

/** Comparison operator */
bool Median::operator!=(const Median &other) const
{
    return not this->operator==(other);
}

/** Completely clear the statistics in this accumulator */
void Median::clear()
{
    minval = 0;
    maxval = 0;
    Accumulator::clear();
}

/** Accumulate the passed value onto the average */
void Median::accumulate(double value)
{
    if (value < minval)
        minval = value;
        
    if (value > maxval)
        maxval = value;

    Accumulator::accumulate(value);
}

/** Return the median value */
double Median::median() const
{
    return (0.5 * maxval) + (0.5 * minval);   // the sum of maxval and minval
                                              // could overflow
}

/** Allow automatic casting to a double to retrieve the average value */
Median::operator double() const
{
    return this->median();
}

/** Return the maximum value */
double Median::max() const
{
    return maxval;
}

/** Return the maximum value */
double Median::maximum() const
{
    return this->max();
}

/** Return the minimum value */
double Median::min() const
{
    return minval;
}

/** Return the minimum value */
double Median::minimum() const
{
    return this->min();
}

const char* Median::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Median>() );
}

/////////
///////// Implementation of RecordValues
/////////

static const RegisterMetaType<RecordValues> r_recval;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const RecordValues &recval)
{
    writeHeader(ds, r_recval, 2);
    
    SharedDataStream sds(ds);
    
    sds << recval.vals
        << static_cast<const Accumulator&>(recval);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, RecordValues &recval)
{
    VersionID v = readHeader(ds, r_recval);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
    
        sds >> recval.vals
            >> static_cast<Accumulator&>(recval);
    }
    else if (v == 1)
    {
        QVector<double> vals;
    
        ds >> vals
           >> static_cast<Accumulator&>(recval);
           
        recval.vals = ChunkedVector<double,2048>::fromVector(vals);
    }
    else
        throw version_error(v, "1", r_recval, CODELOC);
        
    return ds;
}

/** Construct an empty average */
RecordValues::RecordValues() 
             : ConcreteProperty<RecordValues,Accumulator>()
{}

/** Copy constructor */
RecordValues::RecordValues(const RecordValues &other)
             : ConcreteProperty<RecordValues,Accumulator>(other), 
               vals(other.vals)
{}

/** Destructor */
RecordValues::~RecordValues()
{}

/** Copy assignment operator */
RecordValues& RecordValues::operator=(const RecordValues &other)
{
    if (this != &other)
    {
        vals = other.vals;
        Accumulator::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool RecordValues::operator==(const RecordValues &other) const
{
    return vals == other.vals and
           Accumulator::operator==(other);
}

/** Comparison operator */
bool RecordValues::operator!=(const RecordValues &other) const
{
    return not this->operator==(other);
}

/** Completely clear the statistics in this accumulator */
void RecordValues::clear()
{
    vals.clear();
    Accumulator::clear();
}

/** Accumulate the passed value onto the average */
void RecordValues::accumulate(double value)
{
    vals.append(value);

    Accumulator::accumulate(value);
}

/** Return the number of recorded values */
int RecordValues::count() const
{
    return vals.count();
}

/** Return the number of recorded values */
int RecordValues::size() const
{
    return this->count();
}

/** Return the number of recorded values */
int RecordValues::nValues() const
{
    return this->count();
}

/** Return the maximum value */
double RecordValues::max() const
{
    int nvals = this->count();

    if (nvals == 0)
        return 0;
        
    double maxval = -(std::numeric_limits<double>::max());
    
    for (int i=0; i<nvals; ++i)
    {
        if (vals[i] > maxval)
            maxval = vals[i];
    }
    
    return maxval;
}

/** Return the maximum value */
double RecordValues::maximum() const
{
    return this->max();
}

/** Return the minimum value */
double RecordValues::min() const
{
    int nvals = this->count();

    if (nvals == 0)
        return 0;
        
    double minval = std::numeric_limits<double>::max();
    
    for (int i=0; i<nvals; ++i)
    {
        if (vals[i] < minval)
            minval = vals[i];
    }
    
    return minval;
}

/** Return the minimum value */
double RecordValues::minimum() const
{
    return this->min();
}

/** Return the sum of all of the values */
double RecordValues::sum() const
{
    int nvals = vals.count();
    
    double sum = 0;
    
    for (int i=0; i<nvals; ++i)
    {
        sum += vals[i];
    }
    
    return sum;
}

/** Return the sum of the square of all of the values */
double RecordValues::sum2() const
{
    int nvals = vals.count();
    
    double sum2 = 0;
    
    for (int i=0; i<nvals; ++i)
    {
        sum2 += pow_2( vals[i] );
    }
    
    return sum2;
}

/** Return the median value */
double RecordValues::median() const
{
    return (0.5 * min()) + (0.5 * max());   // the sum of maxval and minval
                                            // could overflow
}

/** Return the mean value */
double RecordValues::mean() const
{
    if (this->count() == 0)
        return 0;
        
    else
        return this->sum() / this->count();
}

/** Return the mean of the square values */
double RecordValues::meanOfSquares() const
{
    if (this->count() == 0)
        return 0;
        
    else
        return this->sum2() / this->count();
}

/** Return the standard deviation of the values */
double RecordValues::standardDeviation() const
{
    return std::sqrt( this->meanOfSquares() - pow_2(this->mean()) );
}

/** Return the standard deviation of the values */
double RecordValues::stddev() const
{
    return this->standardDeviation();
}

/** Allow automatic casting to a double to retrieve the mean average value */
RecordValues::operator double() const
{
    return this->mean();
}

/** Return the array of all accumulated values */
QVector<double> RecordValues::values() const
{
    return vals.toVector();
}

const char* RecordValues::typeName()
{
    return QMetaType::typeName( qMetaTypeId<RecordValues>() );
}
