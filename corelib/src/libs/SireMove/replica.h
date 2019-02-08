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

#ifndef SIREMOVE_REPLICA_H
#define SIREMOVE_REPLICA_H

#include "suprasubsystem.h"
#include "ensemble.h"

#include "SireBase/propertymap.h"

#include "SireCAS/symbol.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Replica;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::Replica&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::Replica&);

namespace SireMove
{

class Replicas;
class RepExSubMove;

using SireBase::PropertyName;

using SireCAS::Symbol;

/** This is a replica within a replica exchange simulation
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Replica 
        : public SireBase::ConcreteProperty<Replica,SupraSubSystem>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const Replica&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, Replica&);

friend class Replicas; //so can call protected editing functions
friend class RepExSubMove; //so can call protected editing functions

public:
    Replica();
    Replica(const SupraSubSystem &subsystem);
    
    Replica(const Replica &other);
    
    ~Replica();
    
    Replica& operator=(const Replica &other);
    
    bool operator==(const Replica &other) const;
    bool operator!=(const Replica &other) const;
    
    static const char* typeName();

    const Ensemble& ensemble() const;
    
    const Symbol& energyComponent() const;
    
    const PropertyName& spaceProperty() const;
    
    const Symbol& lambdaComponent() const;
    
    double lambdaValue() const;
    
    SireUnits::Dimension::Temperature temperature() const;
    SireUnits::Dimension::Pressure pressure() const;
    SireUnits::Dimension::Pressure fugacity() const;
    SireUnits::Dimension::MolarEnergy chemicalPotential() const;

    SireUnits::Dimension::Volume volume() const;
    SireUnits::Dimension::MolarEnergy energy();

    bool isConstantEnergy() const;
    bool isConstantTemperature() const;
    bool isConstantPressure() const;
    bool isConstantVolume() const;
    bool isConstantNParticles() const;
    bool isConstantFugacity() const;
    bool isConstantChemicalPotential() const;
    bool isConstantLambda(const Symbol &lam) const;
    
protected:
    void setEnergyComponent(const Symbol &symbol);
    void setSpaceProperty(const PropertyName &spaceproperty);

    void setLambdaComponent(const Symbol &symbol);
    void setLambdaValue(double value);

    void setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void setPressure(const SireUnits::Dimension::Pressure &pressure);
    void setFugacity(const SireUnits::Dimension::Pressure &fugacity);
    void setChemicalPotential(
                     const SireUnits::Dimension::MolarEnergy &chemical_potential);

    void setSubSystem(const System &subsystem);
    void setSubMoves(const Moves &submoves);

    void setSubSystemAndMoves(const SimStore &simstore);

    void setGenerator(const RanGenerator &rangenerator);

    void swapInSystem(const SimStore &simstore, bool swap_monitors=false);
    void swapInMolecules(const SimStore &simstore);

    void _post_unpack();

private:
    void updatedMoves();
    
    /** Identifiers for the list of actions that can be deferred
        until the replica is unpacked */
    enum ReplicaCommand 
         { 
           ENERGY_COMPONENT = 1,  // calls this->setEnergyComponent
           SPACE_PROPERTY   = 2,  // calls this->setSpaceProperty
           LAMBDA_COMPONENT = 3,  // calls this->setLambdaComponent
           LAMBDA_VALUE     = 4,  // calls this->setLambdaValue
           REP_TEMPERATURE  = 5,  // calls this->setTemperature
           REP_PRESSURE     = 6,  // calls this->setPressure
           REP_CHEMPOT      = 7,  // calls this->setChemicalPotential
           REP_FUGACITY     = 8,  // calls this->setFugacity
           SWAP_REP_AND_MON = 9,  // calls this->swapInSystem(replica, true)
           SWAP_REP_ONLY    = 10, // calls this->swapInSystem(replica, false)
           SWAP_MOLECULES   = 11, // calls this->swapInMolecules(replica)
           SET_RANGENERATOR = 12  // calls this->setGenerator
         };

    template<class T>
    void deferCommand(ReplicaCommand command, const T &argument);

    /** The values that need to be set when the replica is unpacked, 
        in the order in which they should be applied */
    QList< QPair<quint32,QVariant> > vars_to_be_set;

    /** The ensemble sampled by moves in this replica */
    Ensemble replica_ensemble;
    
    /** The property used to get the simulation space (for volume)
        for this replica */
    PropertyName space_property;
    
    /** The symbol that represents the component of the energy
        that is evaluate for this replica - it is this energy
        that is put into the replica exchange test */
    Symbol nrg_component;
    
    /** The symbol representing the lambda coordinate for 
        lambda-based Hamiltonian replica exchange. This is null
        if this type of replica-exchange is not being performed */
    Symbol lambda_component;
    
    /** The value of lambda for this replica */
    double lambda_value;
};

}

Q_DECLARE_METATYPE( SireMove::Replica )

SIRE_EXPOSE_CLASS( SireMove::Replica )

SIRE_END_HEADER

#endif
