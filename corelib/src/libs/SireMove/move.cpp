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

#include "move.h"
#include "ensemble.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"

#include "SireMaths/rangenerator.h"

#include "SireSystem/system.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireSystem;
using namespace SireStream;

///////
/////// Implementation of Move
///////

static const RegisterMetaType<Move> r_move( MAGIC_ONLY, "SireMove::Move" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Move &move)
{
    writeHeader(ds, r_move, 2);
    
    SharedDataStream sds(ds);
    
    sds << move.nrgcomponent
        << move.coordsproperty << move.spaceproperty
        << move.map
        << static_cast<const Property&>(move);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Move &move)
{
    VersionID v = readHeader(ds, r_move);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
    
        sds >> move.nrgcomponent
            >> move.coordsproperty >> move.spaceproperty
            >> move.map
            >> static_cast<Property&>(move);
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> move.nrgcomponent
            >> move.coordsproperty >> move.spaceproperty
            >> static_cast<Property&>(move);
            
        move.map = PropertyMap();
    }
    else
        throw version_error(v, "1,2", r_move, CODELOC);
        
    return ds;
}

/** Constructor */
Move::Move(const PropertyMap &map) 
     : Property(), nrgcomponent(ForceFields::totalComponent())
{
    coordsproperty = map["coordinates"];
    spaceproperty = map["space"];
}

/** Copy constructor */
Move::Move(const Move &other) 
     : Property(other), nrgcomponent(other.nrgcomponent),
       coordsproperty(other.coordsproperty), spaceproperty(other.spaceproperty),
       map(other.map)
{}

/** Destructor */
Move::~Move()
{}

/** Copy assignment operator */
Move& Move::operator=(const Move &other)
{
    nrgcomponent = other.nrgcomponent;
    coordsproperty = other.coordsproperty;
    spaceproperty = other.spaceproperty;
    map = other.map;
    
    Property::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool Move::operator==(const Move &other) const
{
    return nrgcomponent == other.nrgcomponent and
           coordsproperty == other.coordsproperty and
           spaceproperty == other.spaceproperty and
           map == other.map;
}

/** Comparison operator */
bool Move::operator!=(const Move &other) const
{
    return not Move::operator==(other);
}

/** Perform a single move on the system 'system' without 
    recording any statistics */
void Move::move(System &system)
{
    this->move(system, 1, false);
}

/** Perform 'n' moves on the system without recording any 
    statistics */
void Move::move(System &system, int nmoves)
{
    this->move(system, nmoves, false);
}

/** Set the energy component that describes the Hamiltonian that this move
    will sample */
void Move::setEnergyComponent(const Symbol &component)
{
    nrgcomponent = component;
}

/** Return the symbol that describes the Hamiltonian that this move 
    will sample */
const Symbol& Move::energyComponent() const
{
    return nrgcomponent;
}

/** Set the property used to find the simulation box for this move. This
    property is only used when the move involves changing the 
    simulation box, or when trying to find the volume using
    the 'volume' function */
void Move::setSpaceProperty(const PropertyName &space_property)
{
    spaceproperty = space_property;
    map.set("space", spaceproperty);
}

/** Return the property used to find the simulation space (box) 
    for this move. This property is only used when the move 
    involves changing the simulation box. */
const PropertyName& Move::spaceProperty() const
{
    return spaceproperty;
}

/** Set the property used to locate the molecule coordinates to be moved */
void Move::setCoordinatesProperty(const PropertyName &coords_property)
{
    coordsproperty = coords_property;
    map.set("coordinates", coordsproperty);
}

/** Return the property used to find the molecule coordinates to be moves */
const PropertyName& Move::coordinatesProperty() const
{
    return coordsproperty;
}

/** Return the property map used to perform this move */
const PropertyMap& Move::propertyMap() const
{
    return map;
}

/** Internal function used to update the value of a property */
void Move::setProperty(const QString &property, const PropertyName &value)
{
    map.set(property, value);
}

/** Return the energy of the system 'system' - this is the energy
    of this configuration that will be seen by this move */
MolarEnergy Move::energy(System &system) const
{
    return system.energy(nrgcomponent);
}

/** Return the volume of the system 'system' - this is the volume
    of this configuration that will be seen by this move */
Volume Move::volume(const System &system) const
{
    return system.property(spaceproperty).asA<Space>().volume();
}

/** Return whether or not this move keeps the total energy constant */
bool Move::isConstantEnergy() const
{
    return this->ensemble().isConstantEnergy();
}

/** Return whether or not this move keeps the temperature constant */
bool Move::isConstantTemperature() const
{
    return this->ensemble().isConstantTemperature();
}

/** Return whether or not this move keeps the volume constant */
bool Move::isConstantVolume() const
{
    return this->ensemble().isConstantVolume();
}

/** Return whether or not this move keeps the pressure constant */
bool Move::isConstantPressure() const
{
    return this->ensemble().isConstantPressure();
}

/** Return whether or not this move keeps the chemical potential constant */
bool Move::isConstantChemicalPotential() const
{
    return this->ensemble().isConstantChemicalPotential();
}

/** Return whether or not this move keeps the fugacity (related to the chemical
    potential) constant */
bool Move::isConstantFugacity() const
{
    return this->ensemble().isConstantFugacity();
}

/** Return the constant temperature that this move samples

    \throw SireError::incompatible_error
*/
Temperature Move::temperature() const
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 does not have a constant temperature as it "
            "samples from the %2")
                .arg( this->what() ).arg(this->ensemble().toString()),
                    CODELOC );
                    
    return this->ensemble().temperature();
}

/** Return the constant pressure that this move samples 

    \throw SireError::incompatible_error
*/
Pressure Move::pressure() const
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 does not have a constant pressure as it "
            "samples from the %2")
                .arg( this->what() ).arg(this->ensemble().toString()),
                    CODELOC );
                    
    return this->ensemble().pressure();
}

/** Return the constant fugacity that this move samples 

    \throw SireError::incompatible_error
*/
Pressure Move::fugacity() const
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 does not have a constant fugacity as it "
            "samples from the %2")
                .arg( this->what() ).arg(this->ensemble().toString()),
                    CODELOC );
                    
    return this->ensemble().fugacity();
}

/** Return the constant chemical potential that this move samples 

    \throw SireError::incompatible_error
*/
MolarEnergy Move::chemicalPotential() const
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 does not have a constant chemical potential as it "
            "samples from the %2")
                .arg( this->what() ).arg(this->ensemble().toString()),
                    CODELOC );
                    
    return this->ensemble().chemicalPotential();
}

/** Return whether or not this move keeps the symbol 'lam' constant */
bool Move::isConstantLambda(const Symbol&) const
{
    return true;
}

/** Set the temperature - this should never be called! */
void Move::_pvt_setTemperature(const Temperature &temperature)
{
    throw SireError::program_bug( QObject::tr(
        "Cannot set the temperature to %1 C in %2, as, despite the move promising "
        "it samples the %3, it doesn't provide a function to set the temperature!")
            .arg(temperature.to(celsius))
            .arg(this->what())
            .arg(this->ensemble().toString()), CODELOC );
}
                        
/** Set the pressure - this should never be called! */
void Move::_pvt_setPressure(const Pressure &pressure)
{
    throw SireError::program_bug( QObject::tr(
        "Cannot set the pressure to %1 atm in %2, as, despite the move promising "
        "it samples the %3, it doesn't provide a function to set the pressure!")
            .arg(pressure.to(atm))
            .arg(this->what())
            .arg(this->ensemble().toString()), CODELOC );
}
                       
/** Set the pressure - this should never be called! */
void Move::_pvt_setFugacity(const Pressure &fugacity)
{
    throw SireError::program_bug( QObject::tr(
        "Cannot set the fugacity to %1 bar in %2, as, despite the move promising "
        "it samples the %3, it doesn't provide a function to set the fugacity!")
            .arg(fugacity.to(bar))
            .arg(this->what())
            .arg(this->ensemble().toString()), CODELOC );
}

/** Set the temperature that this constant temperature move samples 
    to 'temperature'
    
    \throw SireError::incompatible_error
*/
void Move::setTemperature(const Temperature &temperature)
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 cannot not have a constant temperature (%2 C) as it "
            "samples from the %3")
                .arg( this->what() )
                .arg( temperature.to(celsius) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setTemperature(temperature);
}

/** Set the pressure that this constant pressure move samples 
    to 'pressure'
    
    \throw SireError::incompatible_error
*/
void Move::setPressure(const SireUnits::Dimension::Pressure &pressure)
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 cannot not have a constant pressure (%2 atm) as it "
            "samples from the %3")
                .arg( this->what() )
                .arg( pressure.to(atm) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setPressure(pressure);
}

/** Set the chemical potential that this constant chemical potential move samples 
    to 'chemical_potential'
    
    \throw SireError::incompatible_error
*/
void Move::setChemicalPotential(const MolarEnergy &chemical_potential)
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 cannot not have a constant chemical potential "
            "(%2 kcal mol-1) as it "
            "samples from the %3")
                .arg( this->what() )
                .arg(chemical_potential.to(kcal_per_mol) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setFugacity( Ensemble::MuVT(this->temperature(), chemical_potential)
                                        .fugacity() );
}

/** Set the fugacity that this constant fugacity move samples 
    to 'fugacity'
    
    \throw SireError::incompatible_error
*/
void Move::setFugacity(const Pressure &fugacity)
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "The move %1 cannot not have a constant fugacity "
            "(%2 bar) as it "
            "samples from the %3")
                .arg( this->what() )
                .arg(fugacity.to(bar) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setFugacity( fugacity );
}

///////
/////// Implementation of NullMove
///////

static const RegisterMetaType<NullMove> r_nullmove;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullMove &nullmove)
{
    writeHeader(ds, r_nullmove, 1);
    ds << static_cast<const Move&>(nullmove);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullMove &nullmove)
{
    VersionID v = readHeader(ds, r_nullmove);
    
    if (v == 1)
    {
        ds >> static_cast<Move&>(nullmove);
    }
    else
        throw version_error(v, "1", r_nullmove, CODELOC);

    return ds;
}

/** Constructor */
NullMove::NullMove() : ConcreteProperty<NullMove,Move>()
{}

const NullMove& Move::null()
{
    return *(create_shared_null<NullMove>());
}

/** Copy constructor */
NullMove::NullMove(const NullMove &other)
         : ConcreteProperty<NullMove,Move>(other)
{}

/** Destructor */
NullMove::~NullMove()
{}

/** Copy assignment operator */
NullMove& NullMove::operator=(const NullMove &other)
{
    return *this;
}

/** Comparison operator */
bool NullMove::operator==(const NullMove &other) const
{
    return true;
}

/** Comparison operator */
bool NullMove::operator!=(const NullMove &other) const
{
    return false;
}

/** Return a string representation */
QString NullMove::toString() const
{
    return QObject::tr("NullMove");
}

/** There are no statistics to clear */
void NullMove::clearStatistics()
{}

/* The NullMove does not use a random number generator */
void NullMove::setGenerator(const RanGenerator&)
{}

/** NullMove doesn't perform any moves - no matter how hard you try! */
void NullMove::move(System &system, int nmoves, bool record_stats)
{
    return;
}

/** There have been and never will be any NullMove events */
int NullMove::nMoves() const
{
    return 0;
}

/** NullMove doesn't change anything (so must be NVE) */
Ensemble NullMove::ensemble() const
{
    return Ensemble::NVE();
}

const char* NullMove::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullMove>() );
}

NullMove* NullMove::clone() const
{
    return new NullMove(*this);
}

