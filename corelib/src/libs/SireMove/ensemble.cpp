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

#include <cmath>

#include "ensemble.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireMove;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;
using namespace SireMove::detail;

static const RegisterMetaType<Ensemble> r_ensemble;

template<class T>
static QDataStream& operator<<(QDataStream &ds, const CharArray<T> &array)
{
    for (int i=0; i<CharArray<T>::count(); ++i)
    {
        ds << array[i];
    }
    
    return ds;
}

template<class T>
static QDataStream& operator>>(QDataStream &ds, CharArray<T> &array)
{
    for (int i=0; i<CharArray<T>::count(); ++i)
    {
        ds >> array[i];
    }
    
    return ds;
}

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const Ensemble &ensemble)
{
    writeHeader(ds, r_ensemble, 1);
    
    ds << double(ensemble.ensemble_temperature)
       << double(ensemble.ensemble_pressure)
       << double(ensemble.ensemble_fugacity)
       << ensemble.ensemble_state
       << static_cast<const Property&>(ensemble);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, Ensemble &ensemble)
{
    VersionID v = readHeader(ds, r_ensemble);
    
    if (v == 1)
    {
        double temp, press, fug;
        
        ds >> temp >> press >> fug
           >> ensemble.ensemble_state
           >> static_cast<Property&>(ensemble);
        
        ensemble.ensemble_temperature = Temperature(temp);
        ensemble.ensemble_pressure = Pressure(press);
        ensemble.ensemble_fugacity = Pressure(fug);
    }
    else
        throw version_error( v, "1", r_ensemble, CODELOC );
        
    return ds;
}

static const quint8 UNKNOWN = 0; // unknown state

static const quint8 N = 4;       // constant number of particles
static const quint8 Mu = 16;     // constant chemical potential

static const quint8 V = 4;       // constant volume
static const quint8 P = 16;      // constant pressure

static const quint8 E = 4;       // constant total energy
static const quint8 T = 16;      // constant temperature

/** Return the description for the passed state parameters */
static CharArray<quint32> getDescription(quint8 n, quint8 v, quint8 t)
{
    CharArray<quint32> description;
    
    description[0] = n;
    description[1] = v;
    description[2] = t;
    
    return description;
}

/** Return the merged version of the passed two states */
template<class T>
static CharArray<T> merge(const CharArray<T> &des0, 
                          const CharArray<T> &des1)
{
    CharArray<T> merged;
    
    for (int i=0; i<CharArray<T>::count(); ++i)
    {
        merged[i] = qMax( des0[i], des1[i] );
    }
    
    return merged;
}

/** Construct an NVE ensemble */
Ensemble::Ensemble() 
         : ensemble_temperature(0), ensemble_pressure(0),
           ensemble_fugacity(0)
{
    ensemble_state = getDescription( N, V, E );
}

/** Copy constructor */
Ensemble::Ensemble(const Ensemble &other)
         : ensemble_temperature(other.ensemble_temperature),
           ensemble_pressure(other.ensemble_pressure),
           ensemble_fugacity(other.ensemble_fugacity),
           ensemble_state(other.ensemble_state)
{}

/** Destructor */
Ensemble::~Ensemble()
{}

/** Copy assignment operator */
Ensemble& Ensemble::operator=(const Ensemble &other)
{
    if (this != &other)
    {
        ensemble_temperature = other.ensemble_temperature;
        ensemble_pressure = other.ensemble_pressure;
        ensemble_fugacity = other.ensemble_fugacity;
        ensemble_state = other.ensemble_state;
    }
    
    return *this;
}

/** Comparison operator */
bool Ensemble::operator==(const Ensemble &other) const
{
    return ensemble_state == other.ensemble_state and
           ensemble_temperature == other.ensemble_temperature and
           ensemble_pressure == other.ensemble_pressure and
           ensemble_fugacity == other.ensemble_fugacity;
}

/** Comparison operator */
bool Ensemble::operator!=(const Ensemble &other) const
{
    return not this->operator==(other);
}

/** Return whether or not energy is a constant in this ensemble */
bool Ensemble::isConstantEnergy() const
{
    return ensemble_state[2] == E;
}

/** Return whether or not temperature is a constant in this ensemble */
bool Ensemble::isConstantTemperature() const
{
    return ensemble_state[2] == T;
}

/** Return whether or not volume is a constant in this ensemble */
bool Ensemble::isConstantVolume() const
{
    return ensemble_state[1] == V;
}

/** Return whether or not pressure is a constant in this ensemble */
bool Ensemble::isConstantPressure() const
{
    return ensemble_state[1] == P;
}

/** Return whether or not the number of particles is constant in this ensemble */
bool Ensemble::isConstantNParticles() const
{
    return ensemble_state[0] == N;
}

/** Return whether the chemical potential is constant in this ensemble */
bool Ensemble::isConstantFugacity() const
{
    return ensemble_state[0] == Mu;
}

/** Return whether the chemical potential is constant in this ensemble */
bool Ensemble::isConstantChemicalPotential() const
{
    return this->isConstantFugacity();
}

/** Return whether or not this is the NVE ensemble */
bool Ensemble::isNVE() const
{
    return ensemble_state == getDescription( N, V, E );
}

/** Return whether or not this is the NVT ensemble */
bool Ensemble::isNVT() const
{
    return ensemble_state == getDescription( N, V, T );
}

/** Return whether or not this is the NPT ensemble */
bool Ensemble::isNPT() const
{
    return ensemble_state == getDescription( N, P, T );
}

/** Return whether or not this is the MuVT ensemble */
bool Ensemble::isMuVT() const
{
    return ensemble_state == getDescription( Mu, V, T );
}

/** Return whether or not this is the microcanonical (NVE) ensemble */
bool Ensemble::isMicroCanonical() const
{
    return this->isNVE();
}

/** Return whether or not this is the canonical (NVT) ensemble */
bool Ensemble::isCanonical() const
{
    return this->isNVT();
}

/** Return whether or not this is the isothermal-isobaric (NPT) ensemble */
bool Ensemble::isIsothermalIsobaric() const
{
    return this->isNPT();
}

/** Return whether or not this is the grand canonical (MuVT) ensemble */
bool Ensemble::isGrandCanonical() const
{
    return this->isMuVT();
}

/** Return the shorthand string for this ensemble (e.g. NVT) */
QString Ensemble::shortHand() const
{
    QString shorthand(3, ' ');
    
    switch (ensemble_state[0])
    {
        case N:
            shorthand[0] = 'N';
            break;
        case Mu:
            shorthand[0] = 'M';   // #warning Need greek Mu character
            break;
        default:
            shorthand[0] = '?';
            break;
    }

    switch (ensemble_state[1])
    {
        case V:
            shorthand[1] = 'V';
            break;
        case P:
            shorthand[1] = 'P';
            break;
        default:
            shorthand[1] = '?';
            break;
    }
    
    switch (ensemble_state[2])
    {
        case E:
            shorthand[2] = 'E';
            break;
        case T:
            shorthand[2] = 'T';
            break;
        default:
            shorthand[2] = '?';
            break;
    }
    
    return shorthand;
}

/** Return the name of this ensemble (if it has a name) */
QString Ensemble::name() const
{
    if (this->isMicroCanonical())
        return QObject::tr("microcanonical (NVE) ensemble");
        
    else if (this->isCanonical())
        return QObject::tr("canonical (NVT) ensemble");
        
    else if (this->isIsothermalIsobaric())
        return QObject::tr("isothermal-isobaric (NPT) ensemble");
        
    else if (this->isGrandCanonical())
        return QObject::tr("grand canonical (MuVT) ensemble");
        
    else
        return QObject::tr("%1 ensemble").arg(this->shortHand());
}

/** Return a string representation of this ensemble */
QString Ensemble::toString() const
{
    QStringList parts;
    
    if (this->isConstantChemicalPotential())
    {
        parts.append( QObject::tr("chemical potential = %1 kcal mol-1")
                            .arg( this->chemicalPotential().to(kcal_per_mol) ) );
    }
    
    if (this->isConstantPressure())
    {
        parts.append( QObject::tr("pressure = %1 atm")
                            .arg( this->pressure().to(atm) ) );
    }
    
    if (this->isConstantTemperature())
    {
        parts.append( QObject::tr("temperature = %1 C")
                            .arg( this->temperature().to(Celsius()) ) );
    }
    
    if (parts.isEmpty())
        return this->name();
        
    else
        return QString("%1 { %2 }")
                    .arg(this->name(), Sire::toString(parts));
}

/** Return the temperature of this ensemble

    \throw SireError::incompatible_error
*/
Temperature Ensemble::temperature() const
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "The %1 ensemble does not have a constant temperature.")
                .arg(this->shortHand()), CODELOC );

    return ensemble_temperature;
}

/** Return the pressure of this ensemble

    \throw SireError::incompatible_error
*/
Pressure Ensemble::pressure() const
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "The %1 ensemble does not have a constant pressure.")
                .arg(this->shortHand()), CODELOC );

    return ensemble_pressure;
}

/** Return the fugacity of this ensemble 

    \throw SireError::incompatible_error
*/
Pressure Ensemble::fugacity() const
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "The %1 ensemble does not have a constant fugacity.")
                .arg(this->shortHand()), CODELOC );

    return ensemble_fugacity;
}

/** Return the chemical potential of this ensemble

    \throw SireError::incompatible_error
*/
MolarEnergy Ensemble::chemicalPotential() const
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "The %1 ensemble does not have a constant chemical potential.")
                .arg(this->shortHand()), CODELOC );
    
    // mu = mu_0 + RT ln ( f / P_0 )
    //
    //  mu == chemical potential
    //  mu_0 == chemical potential of standard state
    //  f == fugacity
    //  where P_0 is 1 bar (standard state)
    
    // we will actually return mu - mu_0  == RT ln (f / P_0 )
    return MolarEnergy( gasr * ensemble_temperature.value()
                             * std::log( ensemble_fugacity / 1.0*bar ) );
}

/** Merge the two ensembles 'e0' and 'e1' together. This tries to find an ensemble
    that satisfies both, e.g. merging NVE with NVT will give NVT,
    while merging NVE with NPT would give NPT. */
Ensemble Ensemble::merge(const Ensemble &e0, const Ensemble &e1)
{
    Ensemble merged;
    
    merged.ensemble_state = ::merge(e0.ensemble_state, e1.ensemble_state);
                                              
    if (merged.ensemble_state[0] == Mu)
    {
        if (e0.ensemble_state[0] == Mu)
        {
            if (e1.ensemble_state[0] == Mu and
                e0.fugacity() != e1.fugacity())
            {
                merged.ensemble_state[0] = UNKNOWN;
            }
            else
                merged.ensemble_fugacity = e0.fugacity();
        }
        else
            merged.ensemble_fugacity = e1.fugacity();
    }

    if (merged.ensemble_state[1] == P)
    {
        if (e0.ensemble_state[1] == P)
        {
            if (e1.ensemble_state[1] == P and
                e0.pressure() != e1.pressure())
            {
                merged.ensemble_state[1] = UNKNOWN;
            }
            else
                merged.ensemble_pressure = e0.pressure();
        }
        else
            merged.ensemble_pressure = e1.pressure();
    }
    
    if (merged.ensemble_state[2] == T)
    {
        if (e0.ensemble_state[2] == T)
        {
            if (e1.ensemble_state[2] == T and
                e0.temperature() != e1.temperature())
            {
                merged.ensemble_state[2] = UNKNOWN;
            }
            else
                merged.ensemble_temperature = e0.temperature();
        }
        else
            merged.ensemble_temperature = e1.temperature();
    }
    
    return merged;
}

/** Merge this ensemble with 'other'. This tries to find an ensemble
    that satisfies both, e.g. merging NVE with NVT will give NVT,
    while merging NVE with NPT would give NPT. */
Ensemble Ensemble::merge(const Ensemble &other) const
{
    return Ensemble::merge(*this, other);
}

/** Return the NVE ensemble */
Ensemble Ensemble::NVE()
{
    Ensemble nve;
    nve.ensemble_state = getDescription(N, V, E);
    return nve;
}

/** Return the NVT ensemble for the temperature 'temperature' */
Ensemble Ensemble::NVT(const Temperature &temperature)
{
    Ensemble nvt;
    nvt.ensemble_state = getDescription(N, V, T);
    nvt.ensemble_temperature = temperature;
    return nvt;
}

/** Return the NPT ensemble for the temperature 'temperature' and 
    the pressure 'pressure' */
Ensemble Ensemble::NPT(const Temperature &temperature,
                       const Pressure &pressure)
{
    Ensemble npt;
    npt.ensemble_state = getDescription(N, P, T);
    npt.ensemble_temperature = temperature;
    npt.ensemble_pressure = pressure;
    return npt;
}

/** Return the MuVT ensemble for the temperature 'temperature' and
    the fugacity 'fugacity' */
Ensemble Ensemble::MuVT(const Temperature &temperature,
                        const Pressure &fugacity)
{
    Ensemble mupt;
    mupt.ensemble_state = getDescription(Mu, V, T);
    mupt.ensemble_temperature = temperature;
    mupt.ensemble_fugacity = fugacity;
    return mupt;
}
                
/** Return the MuVT ensemble for the temperature 'temperature' and
    the chemical potential 'chemical_potential' */      
Ensemble Ensemble::MuVT(const Temperature &temperature,
                        const MolarEnergy &chemical_potential)
{
    // mu = mu_0 + RT ln ( f / P_0 )
    //
    //  mu == chemical potential
    //  mu_0 == chemical potential of standard state
    //  f == fugacity
    //  where P_0 is 1 bar (standard state)
    
    // we will actually use mu - mu_0  == RT ln (f / P_0 )
    // so chemical potential = mu - mu_0
    return Ensemble::MuVT( temperature,
                           std::exp( chemical_potential / (gasr*temperature) ) * bar );
}
                      
/** Syntactic sugar to return the NVE ensemble */
Ensemble Ensemble::microcanonical()
{
    return Ensemble::NVE();
}

/** Syntactic sugar to return the NVT ensemble */
Ensemble Ensemble::canonical(const Temperature &temperature)
{
    return Ensemble::NVT(temperature);
}

/** Syntactic sugar to return the NPT ensemble */
Ensemble Ensemble::isothermalIsobaric(const Temperature &temperature,
                                      const Pressure &pressure)
{
    return Ensemble::NPT(temperature, pressure);
}

/** Syntactic sugar to return the MuVT ensemble */
Ensemble Ensemble::grandCanonical(const Temperature &temperature,
                                  const Pressure &fugacity)
{
    return Ensemble::MuVT(temperature, fugacity);
}

/** Syntactic sugar to return the MuVT ensemble */
Ensemble Ensemble::grandCanonical(const Temperature &temperature,
                                  const MolarEnergy &chemical_potential)
{
    return Ensemble::MuVT(temperature, chemical_potential);
}

const char* Ensemble::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Ensemble>() );
}

Ensemble* Ensemble::clone() const
{
    return new Ensemble(*this);
}
