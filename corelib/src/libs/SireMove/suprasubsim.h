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

#ifndef SIREMOVE_SUPRASUBSIM_H
#define SIREMOVE_SUPRASUBSIM_H

#include "suprasubsimpacket.h"

#include "SireCluster/node.h"
#include "SireCluster/promise.h"

SIRE_BEGIN_HEADER

namespace SireMove
{

using SireCluster::Node;

/** This class is used to start and manage an active
    sub-simulation of a supra-simulation.
    
    A supra-simulation consists of a collection of
    supra-moves that are applied to a supra-system,
    while a supra-sub-simulation consists of a collection
    of sub-moves that are applied to a sub-system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSubSim
{
public:
    SupraSubSim();
    SupraSubSim(const SupraSubSim &other);
    
    ~SupraSubSim();
    
    SupraSubSim& operator=(const SupraSubSim &other);
    
    bool operator==(const SupraSubSim &other) const;
    bool operator!=(const SupraSubSim &other) const;
    
    static SupraSubSim run(const SupraSubSystem &system, const SupraSubMoves &moves,
                           int nmoves, bool record_stats=true);
                        
    static SupraSubSim run(const SupraSubSystem &system, const SupraSubMove &move,
                           int nmoves, bool record_stats=true);
                        
    static SupraSubSim run(const SupraSubSimPacket &simpacket);
    
    static SupraSubSim run(Node &node,
                           const SupraSubSystem &system, const SupraSubMoves &moves,
                           int nmoves, bool record_stats=true);
                        
    static SupraSubSim run(Node &node,
                           const SupraSubSystem &system, const SupraSubMove &move,
                           int nmoves, bool record_stats=true);

    static SupraSubSim run(Node &node, const SupraSubSimPacket &simpacket);
    
    void abort();
    void stop();
    
    void wait();
    bool wait(int timeout);
    
    bool isRunning();
    
    bool isError();
    void throwError();
    
    bool wasStopped();
    
    bool wasAborted();
    
    bool hasFinished();
    
    float progress();

    SupraSubSimPacket input();
    SupraSubSimPacket interimResult();
    SupraSubSimPacket result();

    SupraSubSystemPtr initialSystem();
    SupraSubMovesPtr initialMoves();
    
    SupraSubSystemPtr interimSystem();
    SupraSubMovesPtr interimMoves();
    
    SupraSubSystemPtr system();
    SupraSubMovesPtr moves();

private:
    SupraSubSim(const SireCluster::Promise &promise);

    /** The promise holding the running sub-simulation */
    SireCluster::Promise sim_promise;
};

}

SIRE_EXPOSE_CLASS( SireMove::SupraSubSim )

SIRE_END_HEADER

#endif
