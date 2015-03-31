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

#ifndef SIREMOVE_SUPRASIM_H
#define SIREMOVE_SUPRASIM_H

#include "suprasimpacket.h"

#include "SireCluster/node.h"
#include "SireCluster/promise.h"

SIRE_BEGIN_HEADER

namespace SireMove
{

using SireCluster::Node;

/** This class is used to start and manage an active
    supra-simulation (a simulation of a SupraSystem).
    
    A supra-simulation consists of a collection of
    supra-moves that are applied to a supra-system.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSim
{
public:
    SupraSim();
    SupraSim(const SupraSim &other);
    
    ~SupraSim();
    
    SupraSim& operator=(const SupraSim &other);
    
    bool operator==(const SupraSim &other) const;
    bool operator!=(const SupraSim &other) const;
    
    static SupraSim run(const SupraSystem &system, const SupraMoves &moves,
                        int nmoves, bool record_stats=true);
                        
    static SupraSim run(const SupraSystem &system, const SupraMove &move,
                        int nmoves, bool record_stats=true);
                        
    static SupraSim run(const SupraSimPacket &simpacket);
    
    static SupraSim run(Node &node,
                        const SupraSystem &system, const SupraMoves &moves,
                        int nmoves, bool record_stats=true);
                        
    static SupraSim run(Node &node,
                        const SupraSystem &system, const SupraMove &move,
                        int nmoves, bool record_stats=true);

    static SupraSim run(Node &node, const SupraSimPacket &simpacket);
    
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

    SupraSimPacket input();
    SupraSimPacket interimResult();
    SupraSimPacket result();

    SupraSystemPtr initialSystem();
    SupraMovesPtr initialMoves();
    
    SupraSystemPtr interimSystem();
    SupraMovesPtr interimMoves();
    
    SupraSystemPtr system();
    SupraMovesPtr moves();

private:
    SupraSim(const SireCluster::Promise &promise);

    /** The promise holding the running supra-simulation */
    SireCluster::Promise sim_promise;
};

}

SIRE_EXPOSE_CLASS( SireMove::SupraSim )

SIRE_END_HEADER

#endif
