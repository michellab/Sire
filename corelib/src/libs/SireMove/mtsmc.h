/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOVE_MTSMC_H
#define SIREMOVE_MTSMC_H

#include "montecarlo.h"
#include "moves.h"

#include "SireSystem/constraints.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class MTSMC;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::MTSMC&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::MTSMC&);

namespace SireMove
{

using SireSystem::Constraint;
using SireSystem::Constraints;

/** This is a multiple-time-step Monte Carlo moves. This uses
    the Metropolis-Hamilton acceptance test to perform
    'nfast' Monte Carlo moves on a fast Hamiltonian,
    and then testing whether the resulting configuration
    is suitable for inclusion in the ensemble generated
    using a slow Hamiltonian.
    
    For an application of this method to enhance sampling
    in QM/MM simulations, see Woods, Manby and Mulholland,
    J. Chem. Phys. 2008.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT MTSMC 
            : public SireBase::ConcreteProperty<MTSMC,MonteCarlo>
{

friend QDataStream& ::operator<<(QDataStream&, const MTSMC&);
friend QDataStream& ::operator>>(QDataStream&, MTSMC&);

public:
    MTSMC(const PropertyMap &map = PropertyMap());

    MTSMC(const Moves &fastmoves, const PropertyMap &map = PropertyMap());
    MTSMC(const Moves &fastmoves, int nfastmoves,
          const PropertyMap &map = PropertyMap());
          
    MTSMC(const Moves &fastmoves, const Symbol &fastcomponent,
          const PropertyMap &map = PropertyMap());
    MTSMC(const Moves &fastmoves, const Symbol &fastcomponent, 
          int nfastmoves, const PropertyMap &map = PropertyMap());
    
    MTSMC(const Moves &fastmoves, const Constraints &slow_constraints,
          const PropertyMap &map = PropertyMap());
    MTSMC(const Moves &fastmoves, const Constraints &slow_constraints,
          int nfastmoves, const PropertyMap &map = PropertyMap());
          
    MTSMC(const Moves &fastmoves, const Symbol &fastcomponent,
          const Constraints &slow_constraints,
          const PropertyMap &map = PropertyMap());
    MTSMC(const Moves &fastmoves, const Symbol &fastcomponent,
          const Constraints &slow_constraints, int nfastmoves,
          const PropertyMap &map = PropertyMap());
    
    MTSMC(const MTSMC &other);
    
    ~MTSMC();
    
    MTSMC& operator=(const MTSMC &other);
    
    static const char* typeName();

    bool operator==(const MTSMC &other) const;
    bool operator!=(const MTSMC &other) const;
    
    QString toString() const;
    
    void setFastMoves(const Moves &fastmoves);
    void setNFastMoves(int nfast);
    
    void setFastEnergyComponent(const Symbol &component);
    void setSlowEnergyComponent(const Symbol &component);

    void addSlowConstraint(const Constraint &constraint);
    void setSlowConstraints(const Constraints &constraints);

    void removeSlowConstraints();
    
    const Constraints& slowConstraints() const;

    const Moves& fastMoves() const;
    int nFastMoves() const;
    
    const Symbol& fastEnergyComponent() const;
    const Symbol& slowEnergyComponent() const;
    
    void clearStatistics();
    
    void setGenerator(const RanGenerator &rangenerator);
    
    void move(System &system, int nmoves, bool record_stats=true);

private:
    /** The collection of fast moves that will be applied to the system */
    MovesPtr fastmoves;
    
    /** The constraints that are applied to the system after each 
        block of fast moves */
    Constraints slow_constraints;
    
    /** The energy component on which the fast moves will operate */
    Symbol fastcomponent;
    
    /** The number of fast moves to apply per slow move */
    quint32 nfastmoves;
};

}

Q_DECLARE_METATYPE( SireMove::MTSMC )

SIRE_EXPOSE_CLASS( SireMove::MTSMC )

SIRE_END_HEADER

#endif
