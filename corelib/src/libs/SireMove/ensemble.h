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

#ifndef SIREMOVE_ENSEMBLE_H
#define SIREMOVE_ENSEMBLE_H

#include "SireBase/property.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Ensemble;
}

QDataStream& operator<<(QDataStream &ds, const SireMove::Ensemble&);
QDataStream& operator>>(QDataStream &ds, SireMove::Ensemble&);

namespace SireMove
{

namespace detail
{

template<class T>
class CharArray
{
public:
    CharArray()
    {
        for (unsigned int i=0; i<sizeof(T); ++i)
        {
            array()[i] = 0;
        }
    }
    
    CharArray(const CharArray &other) : data(other.data)
    {}
    
    ~CharArray()
    {}
    
    CharArray<T>& operator=(const CharArray<T> &other)
    {
        data = other.data;
        return *this;
    }
    
    bool operator==(const CharArray<T> &other) const
    {
        return data == other.data;
    }

    bool operator!=(const CharArray<T> &other) const
    {
        return data != other.data;
    }
    
    static int count()
    {
        return sizeof(T);
    }
    
    quint8& operator[](int i)
    {
        return this->array()[i];
    }
    
    const quint8& operator[](int i) const
    {
        return this->array()[i];
    }

    operator T() const
    {
        return data;
    }

private:
    quint8* array()
    {
        return (quint8*)( &data );
    }
    
    const quint8* array() const
    {
        return (const quint8*)( &data );
    }
    
    T data;
};

}

/** This class describes the ensemble that will be created by 
    a collection of moves (e.g. will the moves sample at constant
    volume or temperature?)
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Ensemble 
           : public SireBase::ConcreteProperty<Ensemble,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const Ensemble&);
friend QDataStream& ::operator>>(QDataStream&, Ensemble&);

public:
    Ensemble();
    
    Ensemble(const Ensemble &other);
    
    ~Ensemble();
    
    Ensemble& operator=(const Ensemble &other);
    
    static const char* typeName();
    
    Ensemble* clone() const;
    
    bool operator==(const Ensemble &other) const;
    bool operator!=(const Ensemble &other) const;
    
    bool isConstantEnergy() const;
    bool isConstantTemperature() const;
    bool isConstantVolume() const;
    bool isConstantPressure() const;

    bool isConstantNParticles() const;
    bool isConstantFugacity() const;
    bool isConstantChemicalPotential() const;
    
    bool isNVE() const;
    bool isNVT() const;
    bool isNPT() const;
    bool isMuVT() const;

    bool isMicroCanonical() const;
    bool isCanonical() const;
    bool isIsothermalIsobaric() const;
    bool isGrandCanonical() const;

    QString toString() const;
    
    QString name() const;
    QString shortHand() const;

    SireUnits::Dimension::Temperature temperature() const;
    SireUnits::Dimension::Pressure pressure() const;
    SireUnits::Dimension::Pressure fugacity() const;
    SireUnits::Dimension::MolarEnergy chemicalPotential() const;

    Ensemble merge(const Ensemble &other) const;
    
    static Ensemble merge(const Ensemble &e0, const Ensemble &e1);

    static Ensemble NVE();
    
    static Ensemble NVT(const SireUnits::Dimension::Temperature &temperature);

    static Ensemble NPT(const SireUnits::Dimension::Temperature &temperature,
                        const SireUnits::Dimension::Pressure &pressure);

    static Ensemble MuVT(const SireUnits::Dimension::Temperature &temperature,
                         const SireUnits::Dimension::Pressure &fugacity);
                          
    static Ensemble MuVT(const SireUnits::Dimension::Temperature &temperature,
                         const SireUnits::Dimension::MolarEnergy &chemical_potential);
                          
    static Ensemble microcanonical();
    
    static Ensemble canonical(const SireUnits::Dimension::Temperature &temperature);
    
    static Ensemble isothermalIsobaric(
                              const SireUnits::Dimension::Temperature &temperature,
                              const SireUnits::Dimension::Pressure &pressure);
    
    static Ensemble grandCanonical(
                              const SireUnits::Dimension::Temperature &temperature,
                              const SireUnits::Dimension::Pressure &fugacity);
    
    static Ensemble grandCanonical(
                            const SireUnits::Dimension::Temperature &temperature,
                            const SireUnits::Dimension::MolarEnergy &chemical_potential);

private:
    /** The fixed temperature of the ensemble (if any) */
    SireUnits::Dimension::Temperature ensemble_temperature;
    
    /** The fixed pressure (if any) */
    SireUnits::Dimension::Pressure ensemble_pressure;
    
    /** The fixed fugacity (if any) */
    SireUnits::Dimension::Pressure ensemble_fugacity;

    /** The state of the ensemble */
    detail::CharArray<quint32> ensemble_state;
};

}

Q_DECLARE_METATYPE( SireMove::Ensemble )

SIRE_EXPOSE_CLASS( SireMove::Ensemble )

SIRE_END_HEADER

#endif
