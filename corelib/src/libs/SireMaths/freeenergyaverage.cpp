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

#include "freeenergyaverage.h"

#include "SireUnits/units.h"

#include "SireMaths/maths.h"
#include "SireMaths/histogram.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireMaths;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of FreeEnergyAverage
////////////

static const RegisterMetaType<FreeEnergyAverage> r_avg;

/** Serialise to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const FreeEnergyAverage &avg)
{
    writeHeader(ds, r_avg, 3);
    
    SharedDataStream sds(ds);
    
    sds << avg.hist
        << avg.is_forwards_free_energy
        << static_cast<const ExpAverage&>(avg);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, FreeEnergyAverage &avg)
{
    VersionID v = readHeader(ds, r_avg);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        sds >> avg.hist >> avg.is_forwards_free_energy >> static_cast<ExpAverage&>(avg);
    }
    else if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> avg.hist >> static_cast<ExpAverage&>(avg);
        
        avg.is_forwards_free_energy = true;
    }
    else if (v == 1)
    {
        ds >> static_cast<ExpAverage&>(avg);
        avg.hist = Histogram();
        avg.is_forwards_free_energy = true;
    }
    else
        throw version_error(v, "1-3", r_avg, CODELOC);
        
    return ds;
}

/** Constructor - this defaults to accumulating the average
    at room temperature (25 C) and collects statistics about the
    free energy using a histogram of bin width 0.5 kcal mol-1 */
FreeEnergyAverage::FreeEnergyAverage()
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * double(25*celsius)) ), hist(0.5)
{}

/** Constructor - this defaults to accumulating the average
    at room temperature (25 C) and collects statistics about the
    free energy using a histogram of bin width 0.5 kcal mol-1, specifying
    whether or not this is a forwards free energy */
FreeEnergyAverage::FreeEnergyAverage(bool forwards)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * double(25*celsius)) ), hist(0.5),
                        is_forwards_free_energy(forwards)
{}

/** Construct an accumulator to accumulate the free energy average
    at the specified temperature, and to collect statistics about
    the free energy using a histogram of bin width 0.5 kcal mol-1 */
FreeEnergyAverage::FreeEnergyAverage(const Temperature &temperature, bool forwards)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * temperature.to(kelvin)) ), hist(0.5),
                        is_forwards_free_energy(forwards)
{}

/** Constructor - this defaults to accumulating the average
    at room temperature (25 C) and collects statistics about the
    free energy using a histogram of the passed bin width. If the binwidth
    is zero, then a histogram of energies is not collected */
FreeEnergyAverage::FreeEnergyAverage(const MolarEnergy &binwidth, bool forwards)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * double(25*celsius)) ), hist(binwidth.value()),
                        is_forwards_free_energy(forwards)
{}

/** Construct an accumulator to accumulate the free energy average
    at the specified temperature, and to collect statistics about
    the free energy using a histogram of passed bin width */
FreeEnergyAverage::FreeEnergyAverage(const Temperature &temperature,
                                     const MolarEnergy &binwidth, bool forwards)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(
                        -1.0 / (k_boltz * temperature.to(kelvin)) ), hist(binwidth.value()),
                        is_forwards_free_energy(forwards)
{}

/** Copy constructor */
FreeEnergyAverage::FreeEnergyAverage(const FreeEnergyAverage &other)
                  : ConcreteProperty<FreeEnergyAverage,ExpAverage>(other),
                    hist(other.hist),
                    is_forwards_free_energy(other.is_forwards_free_energy)
{}

/** Destructor */
FreeEnergyAverage::~FreeEnergyAverage()
{}

/** Copy assignment operator */
FreeEnergyAverage& FreeEnergyAverage::operator=(const FreeEnergyAverage &other)
{
    if (this != &other)
    {
        hist = other.hist;
        is_forwards_free_energy = other.is_forwards_free_energy;
        ExpAverage::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool FreeEnergyAverage::operator==(const FreeEnergyAverage &other) const
{
    return this == &other or
           (hist == other.hist and
            is_forwards_free_energy == other.is_forwards_free_energy and
            ExpAverage::operator==(other));
}

/** Comparison operator */
bool FreeEnergyAverage::operator!=(const FreeEnergyAverage &other) const
{
    return not operator==(other);
}

/** Combine the passed average onto this average */
FreeEnergyAverage& FreeEnergyAverage::operator+=(const FreeEnergyAverage &other)
{
    if (is_forwards_free_energy != other.is_forwards_free_energy)
        throw SireError::incompatible_error( QObject::tr(
                "Cannot combine a 'forwards' free energy with a 'backwards' free energy. "
                "%1 vs. %2").arg(this->toString()).arg(other.toString()),
                    CODELOC );

    ExpAverage::operator+=(other);
    hist += other.hist;
 
    return *this;
}

/** Return the combination of this average plus other */
FreeEnergyAverage FreeEnergyAverage::operator+(const FreeEnergyAverage &other) const
{
    FreeEnergyAverage ret(*this);
    ret += other;
    return ret;
}

/** Return the temperature at which the free energy average
    is being accumulated */
Temperature FreeEnergyAverage::temperature() const
{
    return Temperature( -1.0 / (k_boltz*scaleFactor()) );
}

/** Return the histogram of energies */
const Histogram& FreeEnergyAverage::histogram() const
{
    return hist;
}

QString FreeEnergyAverage::toString() const
{
    return QObject::tr("FreeEnergyAverage( dG = %1 kcal mol-1, average = %2 kcal mol-1 "
                       "stderr = %3 kcal mol-1, "
                       "skew = %4 kcal mol-1, nSamples = %5, "
                       "isForwardsFreeEnergy() = %6 )")
                            .arg(this->average())
                            .arg(histogram().mean())
                            .arg(histogram().standardDeviation())
                            .arg(histogram().skew())
                            .arg(nSamples())
                            .arg(isForwardsFreeEnergy());
}

const char* FreeEnergyAverage::typeName()
{
    return QMetaType::typeName( qMetaTypeId<FreeEnergyAverage>() );
}

/** Clear all data from the accumulator */
void FreeEnergyAverage::clear()
{
    hist = Histogram( hist.binWidth() );
    ExpAverage::clear();
}

/** Accumulate the passed energy difference onto the free energy average */
void FreeEnergyAverage::accumulate(double value)
{
    ExpAverage::accumulate(value);
    hist.accumulate(value);
}

/** Return whether or not this is a forwards free energy */
bool FreeEnergyAverage::isForwardsFreeEnergy() const
{
    return is_forwards_free_energy;
}

/** Return whether or not this is a backwards free energy */
bool FreeEnergyAverage::isBackwardsFreeEnergy() const
{
    return not isForwardsFreeEnergy();
}

/** Return the Taylor series expansion estimate the difference in free energy */
double FreeEnergyAverage::taylorExpansion() const
{
    double dg = hist.mean() - 0.5*k_boltz*temperature() *
                    ( hist.meanOfSquares() - (hist.mean()*hist.mean()) );
    
    if (not is_forwards_free_energy)
        dg *= -1;
    
    return dg;
}

/** Return the average free energy. Note that if this is a backwards free energy,
    then this will return the negative (so that it is easy to combine backwards
    and forwards values) */
double FreeEnergyAverage::fepFreeEnergy() const
{
    double dg = ExpAverage::average();
    
    if (not is_forwards_free_energy)
        dg *= -1;
    
    return dg;
}


/** Return the average free energy. Note that if this is a backwards free energy,
    then this will return the negative (so that it is easy to combine backwards
    and forwards values) */
double FreeEnergyAverage::average() const
{
    return fepFreeEnergy();
}

/** Return the square of the average free energy. Note that if this is a backwards
    free energy, then this will return the negative (so that it is easy to combine
    backwards and forwards values) */
double FreeEnergyAverage::average2() const
{
    double dg2 = ExpAverage::average2();
    
    if (not is_forwards_free_energy)
        dg2 *= -1;
    
    return dg2;
}

FreeEnergyAverage::operator double() const
{
    return average();
}

////////////
//////////// Implementation of BennettsFreeEnergyAverage
////////////

static const RegisterMetaType<BennettsFreeEnergyAverage> r_bennetts;

QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds,
                                         const BennettsFreeEnergyAverage &bennetts)
{
    writeHeader(ds, r_bennetts, 2);
    
    SharedDataStream sds(ds);
    
    sds << bennetts.bennetts_avg << bennetts.bennetts_avg2
        << bennetts.const_offset
        << static_cast<const FreeEnergyAverage&>(bennetts);
    
    return ds;
}

QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, BennettsFreeEnergyAverage &bennetts)
{
    VersionID v = readHeader(ds, r_bennetts);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> bennetts.bennetts_avg >> bennetts.bennetts_avg2
            >> bennetts.const_offset
            >> static_cast<FreeEnergyAverage&>(bennetts);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
        
        MolarEnergy tmp;
        
        sds >> bennetts.bennetts_avg >> tmp
            >> bennetts.bennetts_avg2 >> tmp
            >> static_cast<FreeEnergyAverage&>(bennetts);
        
        bennetts.const_offset = 0*kcal_per_mol;
    }
    else
        throw version_error(v, "1,2", r_bennetts, CODELOC);
    
    return ds;
}

/** Constructor */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage()
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(),
       bennetts_avg(0), bennetts_avg2(0), const_offset(0)
{}

/** Constructor, specifying whether or not this is a forwards free energy */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(0)
{}

/** Construct the average at the specified temperature */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const Temperature &temperature,
                                                     bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature,forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(0)
{}

/** Construct, specifying the value of any constant offset for the ratio (C value)
    and the temperature */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const MolarEnergy &constant,
                                                     const Temperature &temperature,
                                                     bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature,forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(constant)
{}

/** Construct, specifying the constant offset for the ratio (C value) */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const MolarEnergy &constant, bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(constant)
{}

/** Construct at the specified temperature, using a histogram of the specified bin width */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const Temperature &temperature,
                                                     const MolarEnergy &binwidth,
                                                     bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature,binwidth,forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(0)
{}

/** Construct at the specified temperature, using the specificed constant (C value),
    using a histogram of the specified bin width */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const MolarEnergy &constant,
                                                     const Temperature &temperature,
                                                     const MolarEnergy &binwidth,
                                                     bool forwards)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(temperature,binwidth,forwards),
       bennetts_avg(0), bennetts_avg2(0), const_offset(constant)
{}

/** Copy constructor */
BennettsFreeEnergyAverage::BennettsFreeEnergyAverage(const BennettsFreeEnergyAverage &other)
     : ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>(other),
       bennetts_avg(other.bennetts_avg), bennetts_avg2(other.bennetts_avg2),
       const_offset(other.const_offset)
{}

/** Destructor */
BennettsFreeEnergyAverage::~BennettsFreeEnergyAverage()
{}

/** Copy assignment operator */
BennettsFreeEnergyAverage&
BennettsFreeEnergyAverage::operator=(const BennettsFreeEnergyAverage &other)
{
    if (this != &other)
    {
        bennetts_avg = other.bennetts_avg;
        bennetts_avg2 = other.bennetts_avg2;
        const_offset = other.const_offset;
        FreeEnergyAverage::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool BennettsFreeEnergyAverage::operator==(const BennettsFreeEnergyAverage &other) const
{
    return bennetts_avg == other.bennetts_avg and bennetts_avg2 == other.bennetts_avg2 and
           const_offset == other.const_offset and FreeEnergyAverage::operator==(other);
}

/** Comparison operator */
bool BennettsFreeEnergyAverage::operator!=(const BennettsFreeEnergyAverage &other) const
{
    return not operator==(other);
}

/** Return whether or not this is a forwards ratio (the numerator in the expression) */
bool BennettsFreeEnergyAverage::isForwardsRatio() const
{
    return isForwardsFreeEnergy();
}

/** Return whether or not this is a backwards ratio (the denominator in the expression) */
bool BennettsFreeEnergyAverage::isBackwardsRatio() const
{
    return isBackwardsFreeEnergy();
}

/** Return the value of the constant offset to the energy used in the Bennetts average */
MolarEnergy BennettsFreeEnergyAverage::constant() const
{
    return const_offset;
}

/** Self-addition operator */
BennettsFreeEnergyAverage&
BennettsFreeEnergyAverage::operator+=(const BennettsFreeEnergyAverage &other)
{
    if (isForwardsRatio() != other.isForwardsRatio())
        throw SireError::incompatible_error( QObject::tr(
                "Cannot combine a 'forwards' free energy with a 'backwards' free energy. "
                "%1 vs. %2").arg(this->toString()).arg(other.toString()),
                    CODELOC );

    if (this->constant() != other.constant())
        throw SireError::incompatible_error( QObject::tr(
                "Cannot combine Bennetts averages that have different constant offsets! "
                "%1 vs. %2").arg(this->toString()).arg(other.toString()), CODELOC );

    double nsteps = nSamples() + other.nSamples();
        
    double my_ratio = nSamples() / nsteps;
    double other_ratio = other.nSamples() / nsteps;

    FreeEnergyAverage::operator+=(other);
    
    bennetts_avg = bennetts_avg * my_ratio + other.bennetts_avg * other_ratio;
    bennetts_avg2 = bennetts_avg2 * my_ratio + other.bennetts_avg2 * other_ratio;
    
    return *this;
}

/** Addition operator */
BennettsFreeEnergyAverage
BennettsFreeEnergyAverage::operator+(const BennettsFreeEnergyAverage &other) const
{
    BennettsFreeEnergyAverage ret(*this);
    ret += other;
    return ret;
}

const char* BennettsFreeEnergyAverage::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BennettsFreeEnergyAverage>() );
}

QString BennettsFreeEnergyAverage::toString() const
{
    return QObject::tr("BennettsFreeEnergyAverage( dG = %1 kcal mol-1, average = %2 kcal mol-1 "
                       "bennettsRatio() = %3, bennettsStandardError(90) = %4, "
                       "stderr = %5 kcal mol-1, "
                       "skew = %6 kcal mol-1, nSamples = %7, constant() = %8 kcal mol-1 "
                       "isForwardsRatio() = %9 )")
                            .arg(this->fepFreeEnergy())
                            .arg(histogram().mean())
                            .arg(bennettsRatio())
                            .arg(bennettsStandardError(90))
                            .arg(histogram().standardError())
                            .arg(histogram().skew())
                            .arg(nSamples())
                            .arg(constant().to(kcal_per_mol))
                            .arg(isForwardsRatio());
}

/** Clear this accumulator */
void BennettsFreeEnergyAverage::clear()
{
    FreeEnergyAverage::clear();
    bennetts_avg = 0;
    bennetts_avg2 = 0;
}

/** Accumulate the passed value onto the average */
void BennettsFreeEnergyAverage::accumulate(double value)
{
    double nsteps = nSamples() + 1;
        
    double my_ratio = nSamples() / nsteps;
    double other_ratio = 1.0 / nsteps;

    double val;

    // Accumulating the average of    1 / { 1 + e^( beta dE - C) } for forwards,
    //                                1 / { 1 + e^( beta dE + C) } for backwards
    //
    // (note that this->scaleFactor() is -beta)
    
    if (isForwardsRatio())
    {
        val = 1.0 / (1.0 + std::exp(-this->scaleFactor()*(value-const_offset.value())));
    }
    else
    {
        val = 1.0 / (1.0 + std::exp(-this->scaleFactor()*(value+const_offset.value())));
    }

    bennetts_avg = my_ratio * bennetts_avg + other_ratio * val;
    bennetts_avg2 = my_ratio * bennetts_avg2 + other_ratio * (val*val);
    
    FreeEnergyAverage::accumulate(value);
}

/** Return the Bennetts ratio. This is the ensemble average
    of 1 / {1 + exp( beta dE + C ) } if this is a forwards ratio, or
    of 1 / {1 + exp( beta dE - C ) } if this is a backwards ratio */
double BennettsFreeEnergyAverage::bennettsRatio() const
{
    return bennetts_avg;
}

/** Return the standard error on the Bennetts ratio to the passed confidence level */
double BennettsFreeEnergyAverage::bennettsStandardError(double level) const
{
    if (this->nSamples() == 0)
        return 0;

    double stdev = std::sqrt( bennetts_avg2 - pow_2(bennetts_avg) );
    return Histogram::tValue(this->nSamples(),level) * stdev / std::sqrt(this->nSamples());
}
