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

#include "moves.h"
#include "simulation.h"
#include "ensemble.h"

#include "SireMaths/rangenerator.h"

#include "SireUnits/units.h"
#include "SireUnits/temperature.h"
#include "SireUnits/dimensions.h"

#include "SireSystem/system.h"
#include "SireSystem/constraints.h"

#include "SireSystem/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireSystem;
using namespace SireCAS;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireBase;
using namespace SireID;
using namespace SireStream;

///////
/////// Implementation of Moves
///////

static const RegisterMetaType<Moves> r_moves( MAGIC_ONLY, "SireMove::Moves" );
                                                      
/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const Moves &moves)
{
    writeHeader(ds, r_moves, 1);
    
    ds << static_cast<const Property&>(moves);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, Moves &moves)
{
    VersionID v = readHeader(ds, r_moves);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(moves);
    }
    else
        throw version_error(v, "1", r_moves, CODELOC);
        
    return ds;
}

/** Constructor */
Moves::Moves() : Property()
{}

/** Copy constructor */
Moves::Moves(const Moves &other) : Property(other)
{}

/** Destructor */
Moves::~Moves()
{}

/** Return the ith Move object in this set of Moves

    \throw SireError::invalid_index
*/
const Move& Moves::operator[](int i) const
{
    QList<MovePtr> mvs = this->moves();
    
    return mvs.at( Index(i).map(mvs.count()) ).read();
}

/** Return the number of different move types in this set of Moves 
    (this is the maximum number of Move objects that can be iterated over) */
int Moves::nMoveTypes() const
{
    return this->moves().count();
}

/** Return the number of different move types in this set of Moves 
    (this is the maximum number of Move objects that can be iterated over) */
int Moves::count() const
{
    return this->nMoveTypes();
}

/** Return the number of different move types in this set of Moves 
    (this is the maximum number of Move objects that can be iterated over) */
int Moves::size() const
{
    return this->nMoveTypes();
}

/** Return the total number of moves performed using the moves
    in this set */
int Moves::nMoves() const
{
    const QList<MovePtr> mvs = this->moves();
    
    int nmoves = 0;
    
    for (QList<MovePtr>::const_iterator it = mvs.constBegin();
         it != mvs.constEnd();
         ++it)
    {
        nmoves += it->read().nMoves();
    }
    
    return nmoves;
}
    
/** Set the energy component used by these moves - by default 
    this will raise an error as setting the energy component is
    not supported 
    
    \throw SireError::unsupported
*/
void Moves::setEnergyComponent(const Symbol &nrg_component)
{
    throw SireError::unsupported( QObject::tr(
            "The moves object %1 does not support the setting of the "
            "energy component used during sampling (to %2).")
                .arg(this->what())
                .arg(nrg_component.toString()), CODELOC );
}

/** Set the space property used by these moves - by default 
    this will raise an error as setting the space property is
    not supported 
    
    \throw SireError::unsupported
*/
void Moves::setSpaceProperty(const PropertyName &space_property)
{
    throw SireError::unsupported( QObject::tr(
            "The moves object %1 does not support the setting of the "
            "space property used during sampling (to %2).")
                .arg(this->what())
                .arg(space_property.toString()), CODELOC );
}

/** Return the energy of the system that is seen by these moves */
MolarEnergy Moves::energy(System &system) const
{
    return system.energy( this->energyComponent() );
}

/** Return the volume of the system as it is seen by these moves */
Volume Moves::volume(const System &system) const
{
    return system.property( this->spaceProperty() ).asA<Space>().volume();
}

/** Return the ensemble that will be sampled by these moves */
Ensemble Moves::ensemble() const
{
    QList<MovePtr> mvs = this->moves();
    
    if (mvs.isEmpty())
        return Ensemble();
    
    Ensemble merged = mvs.takeFirst()->ensemble();

    for (QList<MovePtr>::const_iterator it = mvs.constBegin();
         it != mvs.constEnd();
         ++it)
    {
        merged = Ensemble::merge(merged, (*it)->ensemble());
    }
    
    return merged;
}

/** Return whether or not these moves keeps the total energy constant */
bool Moves::isConstantEnergy() const
{
    return this->ensemble().isConstantEnergy();
}

/** Return whether or not these moves keep the temperature constant */
bool Moves::isConstantTemperature() const
{
    return this->ensemble().isConstantTemperature();
}

/** Return whether or not these moves keep the volume constant */
bool Moves::isConstantVolume() const
{
    return this->ensemble().isConstantVolume();
}

/** Return whether or not these moves keep the pressure constant */
bool Moves::isConstantPressure() const
{
    return this->ensemble().isConstantPressure();
}

/** Return whether or not these moves keep the chemical potential constant */
bool Moves::isConstantChemicalPotential() const
{
    return this->ensemble().isConstantChemicalPotential();
}

/** Return whether or not these moves keep the fugacity (related to the chemical
    potential) constant */
bool Moves::isConstantFugacity() const
{
    return this->ensemble().isConstantFugacity();
}

/** Return the constant temperature that these moves sample

    \throw SireError::incompatible_error
*/
Temperature Moves::temperature() const
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "These moves do not have a constant temperature as they "
            "sample from the %1")
                .arg(this->ensemble().toString()),
                    CODELOC );
                    
    return this->ensemble().temperature();
}

/** Return the constant pressure that these moves sample 

    \throw SireError::incompatible_error
*/
Pressure Moves::pressure() const
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "Theses moves do not have a constant pressure as they "
            "sample from the %1")
                .arg( this->ensemble().toString() ),
                    CODELOC );
                    
    return this->ensemble().pressure();
}

/** Return the constant fugacity that these moves sample 

    \throw SireError::incompatible_error
*/
Pressure Moves::fugacity() const
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "These moves do not have a constant fugacity as they "
            "sample from the %1")
                .arg( this->ensemble().toString() ),
                    CODELOC );
                    
    return this->ensemble().fugacity();
}

/** Return the constant chemical potential that these moves sample 

    \throw SireError::incompatible_error
*/
MolarEnergy Moves::chemicalPotential() const
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "These moves do not have a constant chemical potential as they "
            "samples from the %1")
                .arg( this->ensemble().toString() ),
                    CODELOC );
                    
    return this->ensemble().chemicalPotential();
}

/** Return whether or not these moves keep the symbol 'lam' constant */
bool Moves::isConstantLambda(const Symbol &lam) const
{
    QList<MovePtr> mvs = this->moves();
    
    for (QList<MovePtr>::const_iterator it = mvs.constBegin();
         it != mvs.constEnd();
         ++it)
    {
        if (not (*it)->isConstantLambda(lam))
            return false;
    }
    
    return true;
}

/** Set the temperature that these constant temperature moves sample 
    to 'temperature'
    
    \throw SireError::incompatible_error
*/
void Moves::setTemperature(const Temperature &temperature)
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "These moves cannot not have a constant temperature (%1 C) as they "
            "sample from the %2")
                .arg( temperature.to(celsius) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setTemperature(temperature);
}

/** Set the pressure that these constant pressure moves sample 
    to 'pressure'
    
    \throw SireError::incompatible_error
*/
void Moves::setPressure(const SireUnits::Dimension::Pressure &pressure)
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "These moves cannot not have a constant pressure (%1 atm) as they "
            "sample from the %2")
                .arg( pressure.to(atm) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setPressure(pressure);
}

/** Set the chemical potential that these constant chemical potential moves sample 
    to 'chemical_potential'
    
    \throw SireError::incompatible_error
*/
void Moves::setChemicalPotential(const MolarEnergy &chemical_potential)
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "These moves cannot not have a constant chemical potential "
            "(%1 kcal mol-1) as they "
            "sample from the %2")
                .arg(chemical_potential.to(kcal_per_mol) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setFugacity( Ensemble::MuVT(this->temperature(), chemical_potential)
                                        .fugacity() );
}

/** Set the fugacity that these constant fugacity moves sample 
    to 'fugacity'
    
    \throw SireError::incompatible_error
*/
void Moves::setFugacity(const Pressure &fugacity)
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "These moves cannot not have a constant fugacity "
            "(%1 bar) as they "
            "samples from the %2")
                .arg(fugacity.to(bar) )
                .arg(this->ensemble().toString()),
                    CODELOC );

    this->_pvt_setFugacity( fugacity );
}

/** Perform a pre-check to ensure that the system is in a sane state
    before the moves are run */
void Moves::preCheck(System &system) const
{
    int ntries = 0;

    while (not system.constraintsSatisfied())
    {
        QStringList unsatisfied_constraints;
    
        //a constraint is violated - get all unsatisfied constraints
        QVector<ConstraintPtr> constraints = system.constraints().constraints();
        
        for (QVector<ConstraintPtr>::const_iterator it = constraints.constBegin();
             it != constraints.constEnd();
             ++it)
        {
            if (not it->read().isSatisfied(system))
            {
                unsatisfied_constraints.append( QObject::tr("%1 : %2")
                        .arg(unsatisfied_constraints.count()+1)
                        .arg(it->read().toString()) );
            }
        }
        
        ++ntries;

        qDebug() << "WARNING - constraints not satisfied\n"
                 << unsatisfied_constraints.join("\n");

        if (ntries > 5)
            throw SireSystem::constraint_error( QObject::tr(
                "Cannot start running the moves as the system contains "
                "unsatisfied constraints. The number of unsatisfied constraints "
                "equals %1.\nHere are the unsatisfied constraints:\n%2")
                    .arg(unsatisfied_constraints.count())
                    .arg(unsatisfied_constraints.join("\n")), CODELOC );

        //apply the constraints now
        system.applyConstraints();
    }
}

/** Perform a post-check to ensure that the system is in a sane state
    after the moves have run */
void Moves::postCheck(System &system) const
{
    if (not system.constraintsSatisfied())
    {
        QStringList unsatisfied_constraints;
    
        //a constraint is violated - get all unsatisfied constraints
        QVector<ConstraintPtr> constraints = system.constraints().constraints();
        
        for (QVector<ConstraintPtr>::const_iterator it = constraints.constBegin();
             it != constraints.constEnd();
             ++it)
        {
            if (not it->read().isSatisfied(system))
            {
                unsatisfied_constraints.append( QObject::tr("%1 : %2")
                        .arg(unsatisfied_constraints.count()+1)
                        .arg(it->read().toString()) );
            }
        }
        
        throw SireSystem::constraint_error( QObject::tr(
                "Performing the moves has resulted in the system containing "
                "unsatisfied constraints. The number of unsatisfied constraints "
                "equals %1.\nHere are the unsatisfied constraints:\n%2\n"
                "Here are the moves that were performed:\n%3")
                    .arg(unsatisfied_constraints.count())
                    .arg(unsatisfied_constraints.join("\n"))
                    .arg(this->toString()), CODELOC );
    }
}

///////
/////// Implementation of SameMoves
///////

static const RegisterMetaType<SameMoves> r_samemoves;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const SameMoves &samemoves)
{
    writeHeader(ds, r_samemoves, 1);
    
    SharedDataStream sds(ds);
    
    sds << samemoves.mv << static_cast<const Moves&>(samemoves);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, SameMoves &samemoves)
{
    VersionID v = readHeader(ds, r_samemoves);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> samemoves.mv >> static_cast<Moves&>(samemoves);
    }
    else
        throw version_error(v, "1", r_samemoves, CODELOC);
        
    return ds;
}

/** Constructor */
SameMoves::SameMoves() : ConcreteProperty<SameMoves,Moves>()
{}

static SharedPolyPointer<SameMoves> shared_null;

const SameMoves& Moves::null()
{
    if (shared_null.constData() == 0)
    {
        QMutexLocker lkr( SireBase::globalLock() );
        
        if (shared_null.constData() == 0)
            shared_null = new SameMoves();
    }

    return *(shared_null.constData());
}

/** Construct to run the move 'move' multiple times */
SameMoves::SameMoves(const Move &move)
          : ConcreteProperty<SameMoves,Moves>(), mv(move)
{}

/** Copy constructor */
SameMoves::SameMoves(const SameMoves &other)
          : ConcreteProperty<SameMoves,Moves>(other), mv(other.mv)
{}

/** Destructor */
SameMoves::~SameMoves()
{}

/** Copy assignment operator */
SameMoves& SameMoves::operator=(const SameMoves &other)
{
    mv = other.mv;
    Moves::operator=(other);
    return *this;
}

/** Comparison operator */
bool SameMoves::operator==(const SameMoves &other) const
{
    return mv == other.mv;
}

/** Comparison operator */
bool SameMoves::operator!=(const SameMoves &other) const
{
    return mv != other.mv;
}

/** Return a string representation of these moves */
QString SameMoves::toString() const
{
    return QObject::tr("SameMoves( %1 )").arg( mv.read().toString() );
}

/** Set the random number generator used at all points in all of the moves */
void SameMoves::setGenerator(const RanGenerator &rangenerator)
{
    mv.edit().setGenerator(rangenerator);
}

/** Apply the move 'nmoves' times to the system 'system', returning
    the result */
System SameMoves::move(const System &system, int nmoves, bool record_stats)
{
    SameMoves old_state(*this);
    System new_system(system);
    
    try
    {
        //ensure that the system is in a sane state before the moves
        this->preCheck(new_system);

        //perform the moves
        mv.edit().move(new_system, nmoves, record_stats);
        
        //ensure that the system has been placed into a sane state
        //after the moves
        this->postCheck(new_system);
    }
    catch(...)
    {
        this->operator=(old_state);
        throw;
    }
    
    return new_system;
}

/** Set the energy component of all of the moves to 'component' */
void SameMoves::setEnergyComponent(const Symbol &component)
{
    if (mv.read().energyComponent() != component)
    {
        mv.edit().setEnergyComponent(component);
    }
}

/** Set the name of the property that all of the moves will use to 
    find the simulation space (simulation box) to 'spaceproperty' */
void SameMoves::setSpaceProperty(const PropertyName &spaceproperty)
{
    if (mv.read().spaceProperty() != spaceproperty)
    {
        mv.edit().setSpaceProperty(spaceproperty);
    }
}

/** Set the temperature for all moves that have a constant temperature
    to 'temperature'. It has already been checked that these moves
    between them sample at constant temperature */
void SameMoves::_pvt_setTemperature(const Temperature &temperature)
{
    if (mv.read().temperature() != temperature)
    {
        mv.edit().setTemperature(temperature);
    }
}

/** Return the space property used to calculate the system volume */
const PropertyName& SameMoves::spaceProperty() const
{
    return mv.read().spaceProperty();
}

/** Return the energy component used to calculate the system energy */
const Symbol& SameMoves::energyComponent() const
{
    return mv.read().energyComponent();
}

/** Set the pressure for all moves that have a constant pressure
    to 'pressure'. It has already been checked that these moves
    between them sample at constant pressure */
void SameMoves::_pvt_setPressure(const Pressure &pressure)
{
    if (mv.read().pressure() != pressure)
    {
        mv.edit().setPressure(pressure);
    }
}

/** Set the fugacity for all moves that have a constant fugacity
    to 'fugacity'. It has already been checked that these moves
    between them sample at constant fugacity */
void SameMoves::_pvt_setFugacity(const Pressure &fugacity)
{
    if (mv.read().fugacity() != fugacity)
    {
        mv.edit().setFugacity(fugacity);
    }
}

/** Clear the move statistics */
void SameMoves::clearStatistics()
{
    mv.edit().clearStatistics();
}

/** Return the moves used by this object */
QList<MovePtr> SameMoves::moves() const
{
    QList<MovePtr> mvs;
    mvs.append(mv);
    
    return mvs;
}

const char* SameMoves::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SameMoves>() );
}

SameMoves* SameMoves::clone() const
{
    return new SameMoves(*this);
}

