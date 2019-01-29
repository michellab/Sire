/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "replica.h"
#include "replicas.h"

#include "SireSystem/system.h"

#include "SireMol/molecules.h"

#include "SireMaths/rangenerator.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireError/errors.h"
#include "SireFF/errors.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<Replica> r_replica;

/** Serialiase to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const Replica &replica)
{
    writeHeader(ds, r_replica, 3);
    
    SharedDataStream sds(ds);
    
    sds << replica.replica_ensemble << replica.space_property
        << replica.nrg_component
        << replica.lambda_component << replica.lambda_value
        << replica.vars_to_be_set
        << static_cast<const SupraSubSystem&>(replica);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, Replica &replica)
{
    VersionID v = readHeader(ds, r_replica);
    
    if (v == 3)
    {
        SharedDataStream sds(ds);
        
        Replica new_replica;
        
        sds >> new_replica.replica_ensemble >> new_replica.space_property
            >> new_replica.nrg_component
            >> new_replica.lambda_component >> new_replica.lambda_value
            >> new_replica.vars_to_be_set
            >> static_cast<SupraSubSystem&>(new_replica);
        
        replica = new_replica;
    }
    else if (v < 3)
    {
        throw SireError::program_bug( QObject::tr(
            "Reading of version %1 of SireMove::Replica should have "
            "been handled by the deprecated classes reader. Something has "
            "gone wrong if they are are being read here...")
                .arg(v), CODELOC );
    }
    else
        throw version_error(v, "3", r_replica, CODELOC);
        
    return ds;
}

static void assertSupportedEnsemble(const Ensemble &replica_ensemble)
{

    if ( not (replica_ensemble.isConstantNParticles() or
              replica_ensemble.isConstantChemicalPotential()) )
        throw SireError::incompatible_error( QObject::tr(
            "Only replica exchange moves involving a constant number of "
            "particles or constant chemical potential are supported. "
            "The %1 is not supported.").arg(replica_ensemble.toString()), 
                CODELOC );

    if ( not (replica_ensemble.isConstantEnergy() or
              replica_ensemble.isConstantTemperature()) )
        throw SireError::incompatible_error( QObject::tr(
            "Only replica exchange moves involving constant energy or "
            "constant temperature ensembles are supported. "
            "The %1 is not supported.").arg(replica_ensemble.toString()), 
                CODELOC );

    if ( not (replica_ensemble.isConstantVolume() or
              replica_ensemble.isConstantPressure()) )
        throw SireError::incompatible_error( QObject::tr(
            "Only replica exchange moves involving constant volume or "
            "constant pressure ensembles are supported. "
            "The %1 is not supported.").arg(replica_ensemble.toString()), 
                CODELOC );
}

/** Internal function called whenever the moves have changed
    in this replica, used to update the replica parameters and 
    ensemble from the moves */
void Replica::updatedMoves()
{
    const Moves &mvs = this->subMoves();

    Symbol mvs_energy = mvs.energyComponent();
    PropertyName mvs_space_property = mvs.spaceProperty();
    Ensemble mvs_ensemble = mvs.ensemble();
    
    ::assertSupportedEnsemble(mvs_ensemble);

    nrg_component = mvs_energy;
    space_property = mvs_space_property;
    replica_ensemble = mvs_ensemble;
    
    if (not lambda_component.isNull())
    {
        if (not mvs.isConstantLambda(lambda_component))
        {
            //these moves don't keep lambda constrained to the same value
            lambda_component = Symbol();
            lambda_value = 0;
        }
    }
}

/** Constructor */
Replica::Replica() : ConcreteProperty<Replica,SupraSubSystem>(), lambda_value(0)
{}

/** Construct from another SubSystem */
Replica::Replica(const SupraSubSystem &subsys)
        : ConcreteProperty<Replica,SupraSubSystem>(subsys), lambda_value(0)
{
    if (subsys.isA<Replica>())
    {
        this->operator=(subsys.asA<Replica>());
    }
    else
    {
        this->updatedMoves();
    }
}

/** Copy constructor */
Replica::Replica(const Replica &other)
        : ConcreteProperty<Replica,SupraSubSystem>(other),
          vars_to_be_set(other.vars_to_be_set),
          replica_ensemble(other.replica_ensemble),
          space_property(other.space_property),
          nrg_component(other.nrg_component),
          lambda_component(other.lambda_component),
          lambda_value(other.lambda_value)
{}

/** Destructor */
Replica::~Replica()
{}

/** Copy assignment operator */
Replica& Replica::operator=(const Replica &other)
{
    if (this != &other)
    {
        vars_to_be_set = other.vars_to_be_set;
        replica_ensemble = other.replica_ensemble;
        space_property = other.space_property;
        nrg_component = other.nrg_component;
        lambda_component = other.lambda_component;
        lambda_value = other.lambda_value;
        
        SupraSubSystem::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool Replica::operator==(const Replica &other) const
{
    return (this == &other) or
           ( replica_ensemble == other.replica_ensemble and
             space_property == other.space_property and
             nrg_component == other.nrg_component and
             lambda_component == other.lambda_component and
             lambda_value == other.lambda_value and
             SupraSubSystem::operator==(other) );
}

/** Comparison operator */
bool Replica::operator!=(const Replica &other) const
{
    return not Replica::operator==(other);
}

/** Return the ensemble defined by the moves of this replica */
const Ensemble& Replica::ensemble() const
{
    return replica_ensemble;
}

/** Return the energy component that describes the Hamiltonian
    that is sampled by this replica */
const Symbol& Replica::energyComponent() const
{
    return nrg_component;
}

/** Return the name of the property containing the simulation box
    that is sampled by this replica - this is used to get the
    volume of the simulation space */
const PropertyName& Replica::spaceProperty() const
{
    return space_property;
}

/** Return the component that can be used to change the Hamiltonian
    for Hamiltonian replica exchange - this is a null symbol if
    this replica is not used in Hamiltonian replica exchange */
const Symbol& Replica::lambdaComponent() const
{
    return lambda_component;
}

/** Return the value of the lambda component - this returns 0 if
    this replica is not suited for Hamiltonian replica exchange */
double Replica::lambdaValue() const
{
    return lambda_value;
}

/** Return the temperature of the replica (if the moves sample 
    a constant temperature ensemble) 
    
    \throw SireError::incompatible_error
*/
Temperature Replica::temperature() const
{
    return replica_ensemble.temperature();
}

/** Return the pressure of the replica (if the moves sample
    a constant pressure ensemble)
    
    \throw SireError::incompatible_error
*/
Pressure Replica::pressure() const
{
    return replica_ensemble.pressure();
}

/** Return the fugacity of the replica (if the moves sample
    a constant chemical potential ensemble) 
    
    \throw SireError::incompatible_error
*/
Pressure Replica::fugacity() const
{
    return replica_ensemble.fugacity();
}

/** Return the chemical potential of the replica (if the moves
    sample a constant chemical potential ensemble)
    
    \throw SireError::incompatible_error
*/
MolarEnergy Replica::chemicalPotential() const
{
    return replica_ensemble.chemicalPotential();
}

/** Return the current volume of the replica's simulation space
    (this could be infinite!) */
Volume Replica::volume() const
{
    const Space &space = this->subSystem().property(space_property).asA<Space>();

    return space.volume();
}

/** Return the total energy of this replica */
MolarEnergy Replica::energy()
{
    System sys = this->subSystem();
    MolarEnergy nrg = sys.energy(nrg_component);
    SupraSubSystem::setSubSystem(sys);
    
    return nrg;
}

/** Return whether or not this is a constant energy replica
    (all moves sample the same total energy) */
bool Replica::isConstantEnergy() const
{
    return replica_ensemble.isConstantEnergy();
}

/** Return whether or not this is a constant temperature replica
    (all moves sample the same temperature) */
bool Replica::isConstantTemperature() const
{
    return replica_ensemble.isConstantTemperature();
}

/** Return whether or not this is a constant pressure replica
    (all moves sample the same pressure) */
bool Replica::isConstantPressure() const
{
    return replica_ensemble.isConstantPressure();
}

/** Return whether or not this is a constant volume replica
    (all moves sample the same volume) */
bool Replica::isConstantVolume() const
{
    return replica_ensemble.isConstantVolume();
}

/** Return whether or not the number of particles is constant
    (all moves keep the same number of particles) */
bool Replica::isConstantNParticles() const
{
    return replica_ensemble.isConstantNParticles();
}

/** Return whether the moves keep the same fugacity */
bool Replica::isConstantFugacity() const
{
    return replica_ensemble.isConstantFugacity();
}

/** Return whether the moves keep the same chemical potential */
bool Replica::isConstantChemicalPotential() const
{
    return replica_ensemble.isConstantChemicalPotential();
}

/** Return whether or not this is a constant lambda replica
    (there is a lambda component, 'lam', and all moves sample the 
     same value of this lambda coordinate) */
bool Replica::isConstantLambda(const Symbol &lam) const
{
    return this->subMoves().isConstantLambda(lam);
}

/** Internal function used to set the sub-system that contains
    the molecular information about this replica. This will
    act to preserve the value of lambda */
void Replica::setSubSystem(const System &subsystem)
{
    if (this->isPacked())
        //need to unpack as we also need to run through 
        //all of the deferred commands
        this->unpack();

    System new_system = subsystem;
    
    if (not lambda_component.isNull())
    {
        new_system.setComponent(lambda_component, lambda_value);
    }
    
    SupraSubSystem::setSubSystem(new_system);
}

/** Internal function used to set the sub-moves that contain
    the move information about this replica. This preserves
    as much ensemble information about the moves as possible,
    e.g. if the original moves are constant temperature and
    the new moves are constant temperature, then the original
    temperature is copied to the new moves */
void Replica::setSubMoves(const Moves &submoves)
{
    if (this->isPacked())
        //need to unpack as we also need to run through
        //all of the deferred commands
        this->unpack();
    
    MovesPtr old_moves = SupraSubSystem::subMoves();
    MovesPtr new_moves = submoves;
    
    try
    {
        if (new_moves->isConstantTemperature() and old_moves->isConstantTemperature())
        {
            if (new_moves->temperature() != old_moves->temperature())
                new_moves.edit().setTemperature( old_moves->temperature() );
        }
    
        if (new_moves->isConstantPressure() and old_moves->isConstantPressure())
        {
            if (new_moves->pressure() != old_moves->pressure())
                new_moves.edit().setPressure( old_moves->pressure() );
        }
    
        if (new_moves->isConstantFugacity() and old_moves->isConstantFugacity())
        {
            if (new_moves->fugacity() != old_moves->fugacity())
                new_moves.edit().setFugacity( old_moves->fugacity() );
        }
    
        if (not nrg_component.isNull())
        {
            if (new_moves->energyComponent() != nrg_component)
                new_moves.edit().setEnergyComponent(nrg_component);
        }
    
        if (not space_property.isNull())
        {
            if (new_moves->spaceProperty() != space_property)
                new_moves.edit().setSpaceProperty(space_property);
        }
    
        SupraSubSystem::setSubMoves(new_moves);
        this->updatedMoves();
    }
    catch(...)
    {
        SupraSubSystem::setSubMoves(old_moves);
        throw;
    }
}

/** Set both the sub-moves and the sub-system that define this replica */
void Replica::setSubSystemAndMoves(const SimStore &simstore)
{
    if (this->isPacked())
        this->unpack();
        
    SimStore unpacked_simstore = simstore;
    
    if (simstore.isPacked())
        unpacked_simstore.unpack();

    //do moves first, as they may change the ensemble
    this->setSubMoves(simstore.moves());
    
    //now do the system
    this->setSubSystem(simstore.system());
}

/** Internal function used to add the command 'command' with 
    argument 'argument' to the list of commands that will be performed
    when the system is next unpacked */
template<class T>
void Replica::deferCommand(ReplicaCommand command, const T &argument)
{
    vars_to_be_set.append( 
        QPair<quint32,QVariant>( quint32(command), QVariant::fromValue<T>(argument) ) );
}

/** Internal function used to set the energy component that represents
    the Hamiltonian sampled by this replica */
void Replica::setEnergyComponent(const Symbol &symbol)
{
    if (this->isPacked())
        this->deferCommand( ENERGY_COMPONENT, symbol );
        
    else
    {
        if (symbol == nrg_component)
            return;
    
        if (not this->subSystem().hasComponent(symbol))
            throw SireFF::missing_component( QObject::tr(
                "Cannot set the energy component for this replica to %1, "
                "as this system (%2) doesn't have such a component. Available energy "
                "components are %3.")
                    .arg(symbol.toString(), subSystem().toString(),
                         Sire::toString(subSystem().componentSymbols())), CODELOC );

        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setEnergyComponent(symbol);
        SupraSubSystem::setSubMoves(mvs);

        nrg_component = symbol;
    }
}

/** Internal function used to set the property used to find the simulation space */
void Replica::setSpaceProperty(const PropertyName &spaceproperty)
{
    if (this->isPacked())
        this->deferCommand( SPACE_PROPERTY, spaceproperty );
        
    else
    {
        if (spaceproperty == space_property)
            return;
    
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setSpaceProperty(spaceproperty);
        SupraSubSystem::setSubMoves(mvs);
        
        space_property = spaceproperty;
    }
}

/** Internal function to set the lambda component used for lambda-based Hamiltonian
    replica exchange */
void Replica::setLambdaComponent(const Symbol &symbol)
{
    if (this->isPacked())
        this->deferCommand(LAMBDA_COMPONENT, symbol);
        
    else
    {
        if (symbol == lambda_component)
            return;
            
        double current_value = 0;
            
        if (not symbol.isNull())
        {
            if (SupraSubSystem::subSystem().hasComponent(symbol))
            {
                System new_system = SupraSubSystem::subSystem();
                current_value = new_system.componentValue(symbol);
            }
            else 
            {
                System new_system = SupraSubSystem::subSystem();
        
                //default always to lambda=0
                new_system.setComponent(symbol, 0);
                current_value = new_system.componentValue(symbol);

                SupraSubSystem::setSubSystem(new_system);
            }
        }
        
        lambda_component = symbol;
        lambda_value = current_value;
    }
}

/** Set the value of the lambda component to 'value' */
void Replica::setLambdaValue(double value)
{
    if (lambda_component.isNull())
        throw SireError::incompatible_error( QObject::tr(
            "You cannot set the value of lambda to %1 as there is no "
            "lambda component.").arg(value), CODELOC );

    if (this->isPacked())
        this->deferCommand(LAMBDA_VALUE, value);
        
    else
    {
        if (value == lambda_value)
            return;
            
        System new_system = SupraSubSystem::subSystem();
        new_system.setComponent( lambda_component, value );
        SupraSubSystem::setSubSystem(new_system);
        
        lambda_value = value;
    }
}

/** Set the temperature of this replica to 'temperature'. This is only possible
    if the moves sample a constant temperature ensemble
    
    \throw SireError::incompatible_error
*/
void Replica::setTemperature(const Temperature &t)
{
    if (not this->isConstantTemperature())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the temperature to %1 K as the replica ensemble (%2) "
            "is not constant temperature.")
                .arg(t.to(kelvin)).arg(replica_ensemble.toString()), CODELOC );
    
    if (this->isPacked())
        this->deferCommand(REP_TEMPERATURE, t.to(kelvin));
        
    else
    {
        if (this->temperature() == t)
            return;
            
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setTemperature(t);
        SupraSubSystem::setSubMoves(mvs);
    }
}

/** Set the pressure of this replica to 'pressure'. This is only possible if
    the moves sample a constant pressure ensemble
    
    \throw SireError::incompatible_error
*/
void Replica::setPressure(const Pressure &p)
{
    if (not this->isConstantPressure())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the pressure to %1 atm as the replica ensemble (%2) "
            "is not constant pressure.")
                .arg(p.to(atm)).arg(replica_ensemble.toString()), CODELOC );

    if (this->isPacked())
        this->deferCommand( REP_PRESSURE, p.to(atm) );

    else
    {
        if (this->pressure() == p)
            return;
        
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setPressure(p);
        SupraSubSystem::setSubMoves(mvs);
    }
}

/** Set the fugacity of this replica to 'fugacity'. This is only possible
    if the moves sample a constant fugacity
    
    \throw SireError::incompatible_error
*/
void Replica::setFugacity(const Pressure &f)
{
    if (not this->isConstantFugacity())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the fugacity to %1 atm as the replica ensemble (%2) "
            "is not constant fugacity.")
                .arg(f.to(atm)).arg(replica_ensemble.toString()), CODELOC );

    if (this->isPacked())
        this->deferCommand( REP_FUGACITY, f.to(atm) );

    else
    {
        if (this->fugacity() == f)
            return;
            
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setFugacity(f);
        SupraSubSystem::setSubMoves(mvs);
    }
}

/** Set the chemical potential of this replica to 'chemical_potential'. 
    This is only possible if the moves sample a constant chemical potential
    
    \throw SireError::incompatible_error
*/
void Replica::setChemicalPotential(const MolarEnergy &c)
{
    if (not this->isConstantChemicalPotential())
        throw SireError::incompatible_error( QObject::tr(
            "Cannot set the chemical potential to %1 kcal mol-1 as the replica "
            "ensemble (%2) is not constant chemical potential.")
                .arg(c.to(kcal_per_mol)).arg(replica_ensemble.toString()), CODELOC );

    if (this->isPacked())
        this->deferCommand( REP_CHEMPOT, c.to(kcal_per_mol) );

    else
    {
        if (this->chemicalPotential() == c)
            return;
            
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setChemicalPotential(c);
        SupraSubSystem::setSubMoves(mvs);
    }
}

/** Set the random number generator used by the replica moves */
void Replica::setGenerator(const RanGenerator &rangenerator)
{
    if (this->isPacked())
        this->deferCommand( SET_RANGENERATOR, rangenerator );
        
    else
    {
        MovesPtr mvs = SupraSubSystem::subMoves();
        mvs.edit().setGenerator(rangenerator);
        SupraSubSystem::setSubMoves(mvs);
    }
}

/** Swap in the sub-system in 'simstore' into this replica. This copies
    the system in 'simstore' into this replica, setting the value
    of lambda of that system to the value for this replica, if 
    necessary. If 'swap_monitors' is true, then the monitors from
    other's system are copied into this replica - otherwise the 
    monitors from the original system are copied into other's system */
void Replica::swapInSystem(const SimStore &simstore, bool swap_monitors)
{
    if (this->isPacked())
    {
        if (swap_monitors)
            this->deferCommand(SWAP_REP_AND_MON, simstore);
        else
            this->deferCommand(SWAP_REP_ONLY, simstore);
    }
    else
    {
        System new_system;

        if (simstore.isPacked())
        {
            SimStore unpacked_simstore = simstore;
            unpacked_simstore.unpack();
        
            new_system = unpacked_simstore.system();
        }
        else
            new_system = simstore.system();
    
        if (not swap_monitors)
            new_system.setMonitors( this->subSystem().monitors() );
            
        this->setSubSystem(new_system);
    }
}

/** Swap in the molecules in 'simstore' into the sub-system that is part
    of this replica. This updates all of the molecules in this system 
    with the status of all of the molecules in 'other'. Note that this 
    won't have any affect unless this replica contains some of the same
    molecules as 'other' */
void Replica::swapInMolecules(const SimStore &simstore)
{
    if (this->isPacked())
        this->deferCommand(SWAP_MOLECULES, simstore);
        
    else
    {
        System old_system;

        if (simstore.isPacked())
        {
            SimStore unpacked_simstore = simstore;
            unpacked_simstore.unpack();
            old_system = unpacked_simstore.system();
        }
        else
            old_system = simstore.system();
        
        System new_system = SupraSubSystem::subSystem();
        new_system.update( old_system.molecules() );
        SupraSubSystem::setSubSystem(new_system);
    }
}

/** Internal function used to extract a value of type 'T' from
    the passed QVariant, throwing an exception if this is not possible
    
    \throw SireError::invalid_cast
*/
template<class T>
static T convert(const QVariant &value)
{
    if (not value.canConvert<T>())
        throw SireError::invalid_cast( QObject::tr(
            "Cannot apply a deferred command as the argument of type %1 "
            "cannot be cast to a value of type %2.")
                .arg( typeid(T).name() )
                .arg( QVariant::typeToName(value.type()) ), CODELOC );
                
    return value.value<T>();
}

/** This function is called just after the system is unpacked.
    This is used to execute commands that have been deferred
    while the system has been packed */
void Replica::_post_unpack()
{
    BOOST_ASSERT( not this->isPacked() );

    //apply the variables in order
    QList< QPair<quint32,QVariant> > to_set = vars_to_be_set;
    
    //clear the list - this is so that we would have the same
    //state if an exception is thrown if the variables had
    //been set directly
    vars_to_be_set = QList< QPair<quint32,QVariant> >();
    
    //apply the actions, in the order that they were requested
    for (QList< QPair<quint32,QVariant> >::const_iterator it = to_set.constBegin();
         it != to_set.constEnd();
         ++it)
    {
        switch (it->first)
        {
            case ENERGY_COMPONENT:
                this->setEnergyComponent( ::convert<Symbol>(it->second) );
                break;
                
            case SPACE_PROPERTY:
                this->setSpaceProperty( ::convert<PropertyName>(it->second) );
                break;
                
            case LAMBDA_COMPONENT:
                this->setLambdaComponent( ::convert<Symbol>(it->second) );
                break;
                
            case LAMBDA_VALUE:
                this->setLambdaValue( ::convert<double>(it->second) );
                break;
                
            case REP_TEMPERATURE:
                this->setTemperature( ::convert<double>(it->second) * kelvin );
                break;
            
            case REP_PRESSURE:
                this->setPressure( ::convert<double>(it->second) * atm );
                break;
                
            case REP_CHEMPOT:
                this->setChemicalPotential( 
                                    ::convert<double>(it->second) * kcal_per_mol  );
                break;
                
            case REP_FUGACITY:
                this->setFugacity( ::convert<double>(it->second) * atm );
                break;
                
            case SWAP_REP_AND_MON:
                this->swapInSystem( ::convert<SimStore>(it->second), true );
                break;
                
            case SWAP_REP_ONLY:
                this->swapInSystem( ::convert<SimStore>(it->second), false );
                break;
                
            case SWAP_MOLECULES:
                this->swapInMolecules( ::convert<SimStore>(it->second) );
                break;
                
            case SET_RANGENERATOR:
                this->setGenerator( ::convert<RanGenerator>(it->second) );
                break;
                
            default:
                throw SireError::unsupported( QObject::tr(
                    "A request was made of an unsuppoted action in Replica. "
                    "The action with ID %1 was requested, but this is not "
                    "supported with this version of Replica.")
                        .arg(it->first), CODELOC );
        }
    }
}

const char* Replica::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Replica>() );
}
