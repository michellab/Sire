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

#ifndef SIREMOVE_REPLICAS_H
#define SIREMOVE_REPLICAS_H

#include "suprasystem.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class Replicas;
}

QDataStream& operator<<(QDataStream&, const SireMove::Replicas&);
QDataStream& operator>>(QDataStream&, SireMove::Replicas&);

namespace SireMaths
{
class RanGenerator;
}

namespace SireMove
{

class Replica;

using SireMaths::RanGenerator;

/** This class is used to hold all of the replicas in 
    a replica exchange simulation
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Replicas 
        : public SireBase::ConcreteProperty<Replicas,SupraSystem>
{

friend QDataStream& ::operator<<(QDataStream&, const Replicas&);
friend QDataStream& ::operator>>(QDataStream&, Replicas&);

public:
    Replicas();
    Replicas(int n);
    Replicas(const System &system, int n=1);
    Replicas(const QVector<System> &systems);
    
    Replicas(const SupraSubSystem &subsystem, int n=1);
    Replicas(const SupraSystem &suprasystem);
    
    Replicas(const Replicas &other);
    
    ~Replicas();
    
    Replicas& operator=(const Replicas &other);
    
    bool operator==(const Replicas &other) const;
    bool operator!=(const Replicas &other) const;
    
    const Replica& operator[](int i) const;
    
    static const char* typeName();

    const Replica& at(int i) const;
    
    int nReplicas() const;
    
    void resetReplicaIDs();
    
    const QVector<quint32>& replicaIDs() const;
    
    void collectSupraStats();
    
    QVector<double> lambdaTrajectory() const;
    
    QList< QVector<double> > lambdaTrajectoryHistory() const;
    
    void setReplicas(const Replicas &replicas);
    
    void setReplica(const Replica &replica);
    void setReplica(int i, const Replica &replica);
    
    void setSubSystem(const System &system);
    void setSubSystem(const SupraSubSystem &subsystem);
    
    void setSubSystem(int i, const System &system);
    void setSubSystem(int i, const SupraSubSystem &subsystem);
    
    void setEnergyComponent(const Symbol &symbol);
    void setEnergyComponent(int i, const Symbol &symbol);
    
    void setSpaceProperty(const PropertyName &spaceproperty);
    void setSpaceProperty(int i, const PropertyName &spaceproperty);
    
    void setLambdaComponent(const Symbol &symbol);
    void setLambdaComponent(int i, const Symbol &symbol);
    
    void setLambdaValue(double value);
    void setLambdaValue(int i, double value);
    
    void setTemperature(const SireUnits::Dimension::Temperature &temperature);
    void setTemperature(int i, const SireUnits::Dimension::Temperature &temperature);
    
    void setPressure(const SireUnits::Dimension::Pressure &pressure);
    void setPressure(int i, const SireUnits::Dimension::Pressure &pressure);
    
    void setFugacity(const SireUnits::Dimension::Pressure &fugacity);
    void setFugacity(int i, const SireUnits::Dimension::Pressure &fugacity);
    
    void setChemicalPotential(
                const SireUnits::Dimension::MolarEnergy &chemical_potential);

    void setChemicalPotential(int i,
                const SireUnits::Dimension::MolarEnergy &chemical_potential);

    void setGenerator(const RanGenerator &rangenerator);
    void setGenerator(int i, const RanGenerator &rangenerator);

    void swapSystems(int i, int j, bool swap_monitors=true);

    void swapMolecules(int i, int j);

protected:
    Replica& _pvt_replica(int i);
    const Replica& _pvt_replica(int i) const;
    
    const Replica& _pvt_constReplica(int i) const;

private:
    /** The index of each of the replicas - this allows the 
        replicas to be tracked as they are swapped around */
    QVector<quint32> replica_ids;
    
    /** The history of lambda values sampled by each replica */
    QList< QVector<double> > replica_history;
};

}

Q_DECLARE_METATYPE( SireMove::Replicas )

SIRE_EXPOSE_CLASS( SireMove::Replicas )

SIRE_END_HEADER

#endif
