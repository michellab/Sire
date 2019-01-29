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

#ifndef SIREMOVE_REPEXMOVE_H
#define SIREMOVE_REPEXMOVE_H

#include "SireMaths/rangenerator.h"

#include "SireMove/supramove.h"
#include "SireMove/suprasubmove.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class RepExSubMove;
class RepExMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::RepExMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::RepExMove&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::RepExSubMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::RepExSubMove&);

namespace SireCluster
{
class Nodes;
}

namespace SireMove
{

class Replicas;
class Replica;

using SireMaths::RanGenerator;

/** This is the sub-move that is applied to each replica in the supra-ensemble

    @author Christopher Woods
*/
class SIREMOVE_EXPORT RepExSubMove
           : public SireBase::ConcreteProperty<RepExSubMove,SupraSubMove>
{

friend QDataStream& ::operator<<(QDataStream&, const RepExSubMove&);
friend QDataStream& ::operator>>(QDataStream&, RepExSubMove&);

public:
    RepExSubMove();
    RepExSubMove(const Replica &replica_a, const Replica &replica_b);
    
    RepExSubMove(const RepExSubMove &other);
    
    ~RepExSubMove();
    
    RepExSubMove& operator=(const RepExSubMove &other);
    
    bool operator==(const RepExSubMove &other) const;
    bool operator!=(const RepExSubMove &other) const;
    
    static const char* typeName();

    QString toString() const;

    SireUnits::Dimension::MolarEnergy energy_i() const;
    SireUnits::Dimension::Volume volume_i() const;
    
    SireUnits::Dimension::MolarEnergy energy_j() const;
    SireUnits::Dimension::Volume volume_j() const;

    void move(SupraSubSystem &system, int n_supra_moves, 
              int n_supra_moves_per_block, bool record_stats);

private:
    void evaluateSwappedState(const Replica &replica);

    template<class T>
    void addPartnerProperty(quint32 property, const T &value);

    /** The volume of the system at the end of the move */
    SireUnits::Dimension::Volume new_volume_i;
    
    /** The energy of the system at the end of the move */
    SireUnits::Dimension::MolarEnergy new_energy_i;
    
    /** The volume of the system at the end of the move in 
        the partner state */
    SireUnits::Dimension::Volume new_volume_j;
    
    /** The energy of the system at the end of the move in
        the new partner state */
    SireUnits::Dimension::MolarEnergy new_energy_j;
    
    enum NewState { LAMBDA_VALUE   = 1,   // the partner has a different lambda value
                    NRG_COMPONENT  = 2,   // the partner samples a different hamiltonian
                    SPACE_PROPERTY = 3    // the partner uses a different space property
                  };
                   
    /** The list of properties of the partner replica that
        this will be swapped with, together with their values */
    QList< QPair<quint32,QVariant> > partner_properties;
    
    /** Whether or not volumes and energies for the new state have been 
        calculated */
    bool have_new_vals;
    
    /** Whether or not the replica move needs the volume of the system */
    bool need_volume;
};

/** This class is used to perform replica exchange moves on a collection
    of Replicas. Each move involves running a block of sampling
    on each of the replicas, and then performing replice exchange swaps
    and tests between pairs.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT RepExMove 
        : public SireBase::ConcreteProperty<RepExMove,SupraMove>
{

friend QDataStream& ::operator<<(QDataStream&, const RepExMove&);
friend QDataStream& ::operator>>(QDataStream&, RepExMove&);

public:
    RepExMove();
    
    RepExMove(const RepExMove &other);
    
    ~RepExMove();
    
    RepExMove& operator=(const RepExMove &other);

    bool operator==(const RepExMove &other) const;
    bool operator!=(const RepExMove &other) const;

    static const char* typeName();
    
    int nAttempted() const;
    int nAccepted() const;
    int nRejected() const;
    
    double acceptanceRatio() const;
    
    void clearStatistics();
    
    void setSwapMonitors(bool swap_monitors);
    
    bool swapMovesDisabled() const;
    void setDisableSwaps(bool disable);
    
    QString toString() const;
    
    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;

    void move(SupraSystem &system, int nmoves, bool record_stats);

private:
    void performMove(SireCluster::Nodes &nodes, Replicas &replicas,
                     bool record_stats);

    bool testPair(const Replica &replica_a, const RepExSubMove &move_a,
                  const Replica &replica_b, const RepExSubMove &move_b) const;

    void testAndSwap(Replicas &replicas, const QVector<RepExSubMove> &submoves,
                     bool even_pairs, bool record_stats);

    /** The random number generator used to accept or reject the moves */
    RanGenerator rangenerator;
    
    /** The number of times a replica exchange move has been accepted */
    quint32 naccept;
    
    /** The number of times a replica exchange move has been rejected */
    quint32 nreject;

    /** Whether or not to swap the system monitors when we swap replicas
         - by default we leave the monitors with the systems */
    bool swap_monitors;
    
    /** Whether or not to disable RETI tests. This is useful when you want
        to just use this to RUN TI on a lot of replicas in parallel */
    bool disable_swaps;
};

}

Q_DECLARE_METATYPE( SireMove::RepExMove )
Q_DECLARE_METATYPE( SireMove::RepExSubMove )

SIRE_EXPOSE_CLASS( SireMove::RepExMove )
SIRE_EXPOSE_CLASS( SireMove::RepExSubMove )

SIRE_END_HEADER

#endif
