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

#ifndef SIREMOVE_SIMULATION_H
#define SIREMOVE_SIMULATION_H

#include "simpacket.h"

#include "SireCluster/node.h"
#include "SireCluster/promise.h"

SIRE_BEGIN_HEADER

namespace SireMove
{

using SireCluster::Node;

/** This class is used start and manage an active
    simulation. A simulation consists of a collection
    of moves that are being applied to a System
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT Simulation
{
public:
    Simulation();
    Simulation(const Simulation &other);

    ~Simulation();

    Simulation& operator=(const Simulation &other);
    
    bool operator==(const Simulation &other) const;
    bool operator!=(const Simulation &other) const;

    static Simulation run( const System &system, const Moves &moves,
                           int nmoves, bool record_stats=true );
                           
    static Simulation run( const System &system, const Moves &moves,
                           int nmoves, int nmoves_per_chunk, 
                           bool record_stats=true );

    static Simulation run( const System &system, const Move &move,
                           int nmoves, bool record_stats=true );
                           
    static Simulation run( const System &system, const Move &move,
                           int nmoves, int nmoves_per_chunk,
                           bool record_stats=true );

    static Simulation run( const SimStore &simstore,
                           int nmoves, bool record_stats=true );
                           
    static Simulation run( const SimStore &simstore,
                           int nmoves, int nmoves_per_chunk, 
                           bool record_stats=true );
                           
    static Simulation run( const SimPacket &simpacket );
                           
    static Simulation run( Node &node,
                           const System &system, const Moves &moves,
                           int nmoves, bool record_stats=true );

    static Simulation run( Node &node,
                           const System &system, const Moves &moves,
                           int nmoves, int nmoves_per_chunk,
                           bool record_stats=true );                           

    static Simulation run( Node &node,
                           const System &system, const Move &move,
                           int nmoves, bool record_stats=true );
                           
    static Simulation run( Node &node,
                           const System &system, const Move &move,
                           int nmoves, int nmoves_per_chunk,
                           bool record_stats=true );
                           
    static Simulation run( Node &node,
                           const SimStore &simstore,
                           int nmoves, bool record_stats=true );

    static Simulation run( Node &node,
                           const SimStore &simstore,
                           int nmoves, int nmoves_per_chunk,
                           bool record_stats=true );                           

    static Simulation run( Node &node, const SimPacket &simpacket );
    
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

    SimPacket input();
    SimPacket interimResult();
    SimPacket result();

    System initialSystem();
    MovesPtr initialMoves();
    
    System interimSystem();
    MovesPtr interimMoves();
    
    System system();
    MovesPtr moves();

private:
    Simulation(const SireCluster::Promise &promise);

    /** The promise holding the running simulation */
    SireCluster::Promise sim_promise;
};

}

SIRE_EXPOSE_CLASS( SireMove::Simulation )

SIRE_END_HEADER

#endif
