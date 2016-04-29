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

#include "moleculardynamics.h"
#include "velocityverlet.h"

#include "SireSystem/system.h"

#include "SireMol/moleculegroup.h"

#include "SireVol/space.h"

#include "SireMaths/rangenerator.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireMol;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

static const RegisterMetaType<MolecularDynamics> r_moldyn;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const MolecularDynamics &moldyn)
{
    writeHeader(ds, r_moldyn, 1);
    
    SharedDataStream sds(ds);
    
    sds << moldyn.intgrator << moldyn.wspace << moldyn.timestep << moldyn.num_moves
        << moldyn.total_time
        << static_cast<const Dynamics&>(moldyn);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        MolecularDynamics &moldyn)
{
    VersionID v = readHeader(ds, r_moldyn);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
    
        sds >> moldyn.intgrator >> moldyn.wspace >> moldyn.timestep 
            >> moldyn.num_moves >> moldyn.total_time
            >> static_cast<Dynamics&>(moldyn);
             
    }
    else
        throw version_error(v, "1", r_moldyn, CODELOC);
        
    return ds;
}

/** Constructor */
MolecularDynamics::MolecularDynamics(const PropertyMap &map)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(map),
                    intgrator( Integrator::null() ), 
                    wspace( IntegratorWorkspace::null() ),
                    timestep(1*femtosecond),
                    num_moves(0), total_time(0)
{
    Dynamics::setEnsemble( Ensemble::NVE() );
}

/** Construct to perform moves on the molecules in the group 'molgroup'. This
    defaults to an all-atom velocity-verlet integrator */
MolecularDynamics::MolecularDynamics(const MoleculeGroup &moleculegroup,
                                     const PropertyMap &map)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(map),
                    intgrator( VelocityVerlet() ), timestep(1*femtosecond), 
                    num_moves(0), total_time(0)
{
    wspace = intgrator.read().createWorkspace(moleculegroup, map);
    Dynamics::setEnsemble( intgrator.read().ensemble() );
}
    
/** Construct a move for the passed molecule group, integrated
    using the supplied integrator */
MolecularDynamics::MolecularDynamics(const MoleculeGroup &moleculegroup, 
                                     const Integrator &integrator,
                                     const PropertyMap &map)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(map),
                    intgrator(integrator), timestep(1*femtosecond), num_moves(0),
                    total_time(0)
{
    wspace = intgrator.read().createWorkspace(moleculegroup, map);
    Dynamics::setEnsemble( intgrator.read().ensemble() );
}

/** Construct a move for the passed molecule group, integrated with 
    the passed timestep */
MolecularDynamics::MolecularDynamics(const MoleculeGroup &molgroup, 
                                     Time t, const PropertyMap &map)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(map),
                    intgrator( VelocityVerlet() ), timestep(t), num_moves(0),
                    total_time(0)
{
    wspace = intgrator.read().createWorkspace(molgroup, map);
    Dynamics::setEnsemble( intgrator.read().ensemble() );
}

/** Construct a move for the passed molecule group, integrated
    using the passed integrator using the passed timestep */
MolecularDynamics::MolecularDynamics(const MoleculeGroup &molgroup,
                                     const Integrator &integrator,
                                     Time t, const PropertyMap &map)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(map),
                    intgrator(integrator), timestep(t), num_moves(0),
                    total_time(0)
{
    wspace = intgrator.read().createWorkspace(molgroup, map);
    Dynamics::setEnsemble( intgrator.read().ensemble() );
}

/** Copy constructor */
MolecularDynamics::MolecularDynamics(const MolecularDynamics &other)
                  : ConcreteProperty<MolecularDynamics,Dynamics>(other),
                    intgrator(other.intgrator), wspace(other.wspace),
                    timestep(other.timestep), num_moves(other.num_moves),
                    total_time(other.total_time)
{}

/** Destructor */
MolecularDynamics::~MolecularDynamics()
{}

/** Copy assignment operator */
MolecularDynamics& MolecularDynamics::operator=(const MolecularDynamics &other)
{
    if (this != &other)
    {
        intgrator = other.intgrator;
        wspace = other.wspace;
        timestep = other.timestep;
        num_moves = other.num_moves;
        total_time = other.total_time;
    
        Dynamics::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool MolecularDynamics::operator==(const MolecularDynamics &other) const
{
    return intgrator == other.intgrator and
           wspace == other.wspace and
           timestep == other.timestep and
           num_moves == other.num_moves and 
           total_time == other.total_time and
           Dynamics::operator==(other);
}

/** Comparison operator */
bool MolecularDynamics::operator!=(const MolecularDynamics &other) const
{
    return not MolecularDynamics::operator==(other);
}

/** Return a string representation of this move */
QString MolecularDynamics::toString() const
{
    return QObject::tr("MolecularDynamics( %1, timeStep() == %2 fs, nMoves() == %3 )")
                .arg(intgrator->toString())
                .arg(timestep.to(femtosecond))
                .arg(num_moves);
}
    
/** Return the number of moves completed using this object */
int MolecularDynamics::nMoves() const
{
    return num_moves;
}
    
/** Return the total amount of time simulated using these moves */
SireUnits::Dimension::Time MolecularDynamics::totalTime() const
{
    return total_time;
}
    
/** Return the molecule group on which this move operates */
const MoleculeGroup& MolecularDynamics::moleculeGroup() const
{
    return wspace.read().moleculeGroup();
}

/** Return the integrator used to advance the coordinates
    from one timestep to the next */
const Integrator& MolecularDynamics::integrator() const
{
    return intgrator.read();
}

/** Set the molecule group containing the molecules to be moved */
void MolecularDynamics::setMoleculeGroup(const MoleculeGroup &new_molgroup)
{
    if (new_molgroup.number() != this->moleculeGroup().number())
    {
        wspace = intgrator.read().createWorkspace(new_molgroup, propertyMap());
    }
}

/** Set the molecule group containing the molecules to be moved */
void MolecularDynamics::setMoleculeGroup(const MoleculeGroup &new_molgroup,
                                         const PropertyMap &map)
{
    if (new_molgroup.number() != this->moleculeGroup().number())
    {
        wspace = intgrator.read().createWorkspace(new_molgroup, map);
    }
    else
    {
        wspace.edit().setPropertyMap(map);
    }
}

/** Set the integrator to be used to advance the coordinates from 
    one timestep to the next. */
void MolecularDynamics::setIntegrator(const Integrator &integrator)
{
    if (intgrator != integrator)
    {
        intgrator = integrator;
        wspace = intgrator.read().createWorkspace( moleculeGroup(), propertyMap() );
        Dynamics::setEnsemble( intgrator.read().ensemble() );
    }
}

/** Set the property used to find the coordinates of the molecules */
void MolecularDynamics::setCoordinatesProperty(const PropertyName &value)
{
    wspace.edit().setCoordinatesProperty(value);
    Move::setCoordinatesProperty(value);
}

/** Set the property used to find the system space */
void MolecularDynamics::setSpaceProperty(const PropertyName &value)
{
    wspace.edit().setSpaceProperty(value);
    Move::setSpaceProperty(value);
}

/** Set the property used to find the molecular velocities */
void MolecularDynamics::setVelocitiesProperty(const PropertyName &value)
{
    wspace.edit().setVelocitiesProperty(value);
    Move::setProperty("velocity", value);
}

/** Set the property used to find the molecular masses */
void MolecularDynamics::setMassesProperty(const PropertyName &value)
{
    wspace.edit().setMassesProperty(value);
    Move::setProperty("mass", value);
}

/** Set the property used to find the elements of the atoms */
void MolecularDynamics::setElementsProperty(const PropertyName &value)
{
    wspace.edit().setElementsProperty(value);
    Move::setProperty("element", value);
}

/** Set the property used to find the generator used to
    generate velocities when they are missing */
void MolecularDynamics::setVelocityGeneratorProperty(const PropertyName &value)
{
    wspace.edit().setVelocityGeneratorProperty(value);
    Move::setProperty("velocity generator", value);
}

/** Return the property used to find the molecular coordinates */
PropertyName MolecularDynamics::coordinatesProperty() const
{
    return propertyMap()["coordinates"];
}

/** Return the property used to find the system space */
PropertyName MolecularDynamics::spaceProperty() const
{
    return propertyMap()["space"];
}

/** Return the property used to find the molecular velocities */
PropertyName MolecularDynamics::velocitiesProperty() const
{
    return propertyMap()["velocity"];
}

/** Return the property used to find the molecular masses */
PropertyName MolecularDynamics::massesProperty() const
{
    return propertyMap()["mass"];
}

/** Return the property used to find the atomic elements */
PropertyName MolecularDynamics::elementsProperty() const
{
    return propertyMap()["element"];
}

/** Return the property used to find the generator for 
    missing velocities */
PropertyName MolecularDynamics::velocityGeneratorProperty() const
{
    return propertyMap()["velocity generator"];
}

/** Return the timestep for the integration */
Time MolecularDynamics::timeStep() const
{
    return timestep;
}

/** Set the timestep for the dynamics integration */
void MolecularDynamics::setTimeStep(const Time &t)
{
    timestep = t;
}

/** Return the kinetic energy of the system at the last move. */
MolarEnergy MolecularDynamics::kineticEnergy() const
{
    return wspace.read().kineticEnergy();
}

/** Return the temperature of the system at the last move */
Temperature MolecularDynamics::temperature() const
{
  SireUnits::Dimension::MolarEnergy ekin = MolecularDynamics::kineticEnergy();

  // NOTE THAT THIS ONLY WORKS FOR 3D SPACE WHEN THERE IS ONE ATOM PER MOLECULE 
  // AND NO CONSTRAINTS...IN OTHER WORDS..NEED TO FIX THIS
  //int ndofs = 3 * wspace.read().nMolecules();
  int ndofs = 3 * 256;

  SireUnits::Dimension::Temperature temp = ( ( 2 * ekin.value() ) / ( ndofs * k_boltz ) ) * kelvin ;
  
  return temp;

}

/** Completely clear any move statistics - this clears all existing
    velocities */
void MolecularDynamics::clearStatistics()
{
    wspace = intgrator.read().createWorkspace( moleculeGroup(), propertyMap() );
    num_moves = 0;
    total_time = Time(0);
}

/** Set the random number generator used by this move
    (this move may be completely deterministic, so may not
     use a generator) */
void MolecularDynamics::setGenerator(const RanGenerator &generator)
{
    wspace.edit().setGenerator(generator);
}

/** Regenerate all of the velocities using the passed velocity generator */
void MolecularDynamics::regenerateVelocities(const System &system,
                                             const VelocityGenerator &generator)
{
    wspace.edit().setSystem(system);
    wspace.edit().regenerateVelocities(generator);
}

/** Perform this move on the System 'system' - perform the move
    'nmoves' times, optionally recording simulation statistics
    if 'record_stats' is true */
void MolecularDynamics::move(System &system, int nmoves, bool record_stats)
{
    if (nmoves <= 0)
        return;

    //first, calculate the energy of the system - this is to prevent a
    //weird bug that causes systems to blow up if the energy has not been
    //calculated - THIS NEEDS FIXING!
    system.energy();

    MolecularDynamics old_state(*this);
        
    System old_system_state(system);
    
    try
    {
        wspace.edit().setSystem(system);

        intgrator.edit().integrate(wspace.edit(), this->energyComponent(),
                                   timestep, nmoves, record_stats);

        system = wspace.read().system();

        num_moves += nmoves;
        total_time += nmoves*timestep;
    }
    catch(...)
    {
        system = old_system_state;
        this->operator=(old_state);

        throw;
    }
}

const char* MolecularDynamics::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MolecularDynamics>() );
}

MolecularDynamics* MolecularDynamics::clone() const
{
    return new MolecularDynamics(*this);
}
