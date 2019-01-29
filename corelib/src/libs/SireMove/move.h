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

#ifndef SIREMOVE_MOVE_H
#define SIREMOVE_MOVE_H

#include "SireCAS/symbol.h"

#include "SireBase/property.h"
#include "SireBase/propertymap.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Move;
class NullMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::Move&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::Move&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullMove&);

namespace SireSystem
{
class System;
}

namespace SireMol
{
class PartialMolecule;
class Molecules;
}

namespace SireMaths
{
class RanGenerator;
}

namespace SireMove
{

class Ensemble;

using SireCAS::Symbol;

using SireMaths::RanGenerator;

using SireSystem::System;

using SireBase::PropertyName;
using SireBase::PropertyMap;

using SireMol::PartialMolecule;
using SireMol::Molecules;

/** This is the base class of all of the move classes

    @author Christopher Woods
*/
class SIREMOVE_EXPORT Move : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const Move&);
friend QDataStream& ::operator>>(QDataStream&, Move&);

public:
    Move(const PropertyMap &map = PropertyMap());
    
    Move(const Move &other);
    
    virtual ~Move();
    
    static const char* typeName()
    {
        return "SireMove::Move";
    }
    
    virtual Move* clone() const=0;
    
    virtual void clearStatistics()=0;
    
    virtual QString toString() const=0;

    virtual int nMoves() const=0;
    
    virtual void move(System &system, int nmoves, bool record_stats)=0;

    void move(System &system);
    void move(System &system, int nmoves);

    virtual SireUnits::Dimension::MolarEnergy energy(System &system) const;
    virtual SireUnits::Dimension::Volume volume(const System &system) const;

    virtual void setGenerator(const RanGenerator &rangenerator)=0;

    const Symbol& energyComponent() const;
    virtual void setEnergyComponent(const Symbol &component);

    const PropertyName& coordinatesProperty() const;
    const PropertyName& spaceProperty() const;

    virtual void setCoordinatesProperty(const PropertyName &coords_property);
    virtual void setSpaceProperty(const PropertyName &space_property);

    const PropertyMap& propertyMap() const;

    virtual Ensemble ensemble() const=0;

    bool isConstantEnergy() const;
    bool isConstantTemperature() const;
    bool isConstantVolume() const;
    bool isConstantPressure() const;
    bool isConstantChemicalPotential() const;
    bool isConstantFugacity() const;
    
    virtual bool isConstantLambda(const Symbol &lam) const;

    SireUnits::Dimension::Temperature temperature() const;
    SireUnits::Dimension::Pressure pressure() const;
    SireUnits::Dimension::Pressure fugacity() const;
    SireUnits::Dimension::MolarEnergy chemicalPotential() const;

    void setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void setPressure(const SireUnits::Dimension::Pressure &pressure);
    void setChemicalPotential(
                    const SireUnits::Dimension::MolarEnergy &chemical_potential);
    void setFugacity(const SireUnits::Dimension::Pressure &fugacity);

    static const NullMove& null();

protected:
    Move& operator=(const Move &other);

    bool operator==(const Move &other) const;
    bool operator!=(const Move &other) const;
    
    void setProperty(const QString &property, const PropertyName &value);
    
    virtual void _pvt_setTemperature(
                            const SireUnits::Dimension::Temperature &temperature);
                            
    virtual void _pvt_setPressure(
                            const SireUnits::Dimension::Pressure &pressure);
                            
    virtual void _pvt_setFugacity(
                            const SireUnits::Dimension::Pressure &fugacity);

private:
    /** The component of the energy that describes the Hamiltonian
        that this move samples */
    Symbol nrgcomponent;
    
    /** The name of the property that contains the Molecule coordinates
        property. This move will only affect the coordinates that are
        contained in this property */
    PropertyName coordsproperty;
    
    /** The name of the property that contains the System space property.
        This is necessary as we may have to map the molecules back into
        the space at the end of the move */
    PropertyName spaceproperty;

    /** The property map used for this move */
    PropertyMap map;
};

/** This is a null move - it doesn't change the system at all! */
class SIREMOVE_EXPORT NullMove : public SireBase::ConcreteProperty<NullMove,Move>
{
public:
    NullMove();
    NullMove(const NullMove &other);
    
    ~NullMove();
    
    NullMove& operator=(const NullMove &other);
    
    bool operator==(const NullMove &other) const;
    bool operator!=(const NullMove &other) const;
    
    static const char* typeName();
    
    NullMove* clone() const;
    
    QString toString() const;
    
    int nMoves() const;
    
    void clearStatistics();
    
    void setGenerator(const RanGenerator &rangenerator);
    
    void move(System &system, int nmoves, bool record_stats);
    
    Ensemble ensemble() const;
};

typedef SireBase::PropPtr<Move> MovePtr;

}

Q_DECLARE_METATYPE( SireMove::NullMove )

SIRE_EXPOSE_CLASS( SireMove::Move )
SIRE_EXPOSE_CLASS( SireMove::NullMove )

SIRE_EXPOSE_PROPERTY( SireMove::MovePtr, SireMove::Move )

SIRE_END_HEADER

#endif
