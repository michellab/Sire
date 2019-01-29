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

#ifndef SIREMATHS_FREEENERGYAVERAGE_H
#define SIREMATHS_FREEENERGYAVERAGE_H

#include "accumulator.h"
#include "histogram.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/temperature.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class FreeEnergyAverage;
class BennettsFreeEnergyAverage;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::FreeEnergyAverage&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::FreeEnergyAverage&);

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::BennettsFreeEnergyAverage&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::BennettsFreeEnergyAverage&);

namespace SireMaths
{

/** This class is used to accumulate the free energy
    average. As well as calculating the average, it also
    records a histogram of values that can be used for
    error analysis
    
    @author Christopher Woods
*/
class SIREMATHS_EXPORT FreeEnergyAverage
           : public SireBase::ConcreteProperty<FreeEnergyAverage,ExpAverage>
{

friend QDataStream& ::operator<<(QDataStream&, const FreeEnergyAverage&);
friend QDataStream& ::operator>>(QDataStream&, FreeEnergyAverage&);

public:
    FreeEnergyAverage();
    FreeEnergyAverage(bool forwards_free_energy);
    FreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                      bool forwards_free_energy=true);
    
    FreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &binwidth,
                      bool forwards_free_energy=true);
    
    FreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                      const SireUnits::Dimension::MolarEnergy &binwidth,
                      bool forwards_free_energy=true);
    
    FreeEnergyAverage(const FreeEnergyAverage &other);
    
    ~FreeEnergyAverage();
    
    FreeEnergyAverage& operator=(const FreeEnergyAverage &other);
    
    bool operator==(const FreeEnergyAverage &other) const;
    bool operator!=(const FreeEnergyAverage &other) const;
    
    FreeEnergyAverage operator+(const FreeEnergyAverage &other) const;
    
    FreeEnergyAverage& operator+=(const FreeEnergyAverage &other);
    
    static const char* typeName();

    SireUnits::Dimension::Temperature temperature() const;

    QString toString() const;

    const Histogram& histogram() const;

    bool isForwardsFreeEnergy() const;
    bool isBackwardsFreeEnergy() const;

    void clear();
    
    void accumulate(double value);

    double average() const;
    double average2() const;

    operator double() const;

    double fepFreeEnergy() const;
    double taylorExpansion() const;

private:
    /** A histogram of all of the energy values. This allows the free energy
        average to be monitored and the error, and corrected values to be 
        calculated */
    Histogram hist;
    
    /** Whether or not this is a forwards free energy (from a low lambda to high lambda
        value) or a backwards free energy (from a high lambda to low lambda). Note that
        the negative of the free energy is returned if this is a backwards free energy */
    bool is_forwards_free_energy;
};

/** This class is used to accumulate the free energy average, using
    both FEP and Bennett's acceptance ratio method. Bennett's method
    calculates a free energy difference between A and B using 
    simulations at both ensembles A and B. The forwards average
    < f( beta(U_B-U_A) )>_A divided by the backwards average
    < f( beta(U_A-U_B) )>_B is equal to exp(-beta deltaG).
    
    In this case, we don't use the "C" energy offset value, as this
    is optimally equal to deltaG, which we don't know, and also because
    in windowed calculations, abs(deltaG) is < 5 kcal mol-1.
    
    Also note that we use the function f(x) = 1 / (1 + exp(x)) as
    this is also supposed to be the most efficient. 
 
    @author Christopher Woods
*/
class SIREMATHS_EXPORT BennettsFreeEnergyAverage
       : public SireBase::ConcreteProperty<BennettsFreeEnergyAverage,FreeEnergyAverage>
{

friend QDataStream& ::operator<<(QDataStream&, const BennettsFreeEnergyAverage&);
friend QDataStream& ::operator>>(QDataStream&, BennettsFreeEnergyAverage&);

public:
    BennettsFreeEnergyAverage();
    BennettsFreeEnergyAverage(bool forwards_free_energy);

    BennettsFreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                              bool forwards_free_energy=true);

    BennettsFreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &constant,
                              const SireUnits::Dimension::Temperature &temperature,
                              bool forwards_free_energy=true);
    
    BennettsFreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &constant,
                              bool forwards_free_energy=true);
    
    BennettsFreeEnergyAverage(const SireUnits::Dimension::Temperature &temperature,
                              const SireUnits::Dimension::MolarEnergy &binwidth,
                              bool forwards_free_energy=true);
    
    BennettsFreeEnergyAverage(const SireUnits::Dimension::MolarEnergy &constant,
                              const SireUnits::Dimension::Temperature &temperature,
                              const SireUnits::Dimension::MolarEnergy &binwidth,
                              bool forwards_free_energy=true);
    
    BennettsFreeEnergyAverage(const BennettsFreeEnergyAverage &other);
    
    ~BennettsFreeEnergyAverage();
    
    BennettsFreeEnergyAverage& operator=(const BennettsFreeEnergyAverage &other);
    
    bool operator==(const BennettsFreeEnergyAverage &other) const;
    bool operator!=(const BennettsFreeEnergyAverage &other) const;
    
    BennettsFreeEnergyAverage operator+(const BennettsFreeEnergyAverage &other) const;
    
    BennettsFreeEnergyAverage& operator+=(const BennettsFreeEnergyAverage &other);
    
    static const char* typeName();

    QString toString() const;

    void clear();
    
    void accumulate(double value);

    double bennettsRatio() const;
    double bennettsStandardError(double level) const;

    bool isForwardsRatio() const;
    bool isBackwardsRatio() const;

    SireUnits::Dimension::MolarEnergy constant() const;

private:
    /** The average bennetts ratio */
    double bennetts_avg;
    
    /** The average of the squared bennetts ratio */
    double bennetts_avg2;
    
    /** The energy offset for each ratio */
    SireUnits::Dimension::MolarEnergy const_offset;
};

}

Q_DECLARE_METATYPE( SireMaths::FreeEnergyAverage )
Q_DECLARE_METATYPE( SireMaths::BennettsFreeEnergyAverage )

SIRE_EXPOSE_CLASS( SireMaths::FreeEnergyAverage )
SIRE_EXPOSE_CLASS( SireMaths::BennettsFreeEnergyAverage )

SIRE_END_HEADER

#endif
