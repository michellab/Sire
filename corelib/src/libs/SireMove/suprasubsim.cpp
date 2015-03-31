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

#include "suprasubsim.h"
#include "suprasubmoves.h"

#include "SireCluster/nodes.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireCluster;

/** Null constructor */
SupraSubSim::SupraSubSim()
{}

/** Internal constructor used to connect to a running simulation */
SupraSubSim::SupraSubSim(const SireCluster::Promise &promise)
            : sim_promise(promise)
{}

/** Copy constructor */
SupraSubSim::SupraSubSim(const SupraSubSim &other)
            : sim_promise(other.sim_promise)
{}

/** Destructor */
SupraSubSim::~SupraSubSim()
{}

/** Copy assignment operator */
SupraSubSim& SupraSubSim::operator=(const SupraSubSim &other)
{
    sim_promise = other.sim_promise;
    return *this;
}

/** Comparison operator */
bool SupraSubSim::operator==(const SupraSubSim &other) const
{
    return sim_promise == other.sim_promise;
}

/** Comparison operator */
bool SupraSubSim::operator!=(const SupraSubSim &other) const
{
    return sim_promise != other.sim_promise;
}

/** Run the sub-system simulation described in 'simpacket' on the node 'node' 
    and return a handle to the running simulation */
SupraSubSim SupraSubSim::run(Node &node, const SupraSubSimPacket &simpacket)
{
    return SupraSubSim( node.startJob(simpacket) );
}
                    
/** Run the sub-system simulation described
    in 'simpacket' in the current thread */
SupraSubSim SupraSubSim::run(const SupraSubSimPacket &simpacket)
{
    Nodes nodes;
    
    ThisThread this_thread = nodes.borrowThisThread();
    
    if (nodes.isEmpty())
        throw SireError::unavailable_resource( QObject::tr(
            "This thread is unavailable for running a simulation. It is already "
            "busy doing something else!"), CODELOC );
            
    Node node = nodes.getNode();
    
    SupraSubSim sim = SupraSubSim::run(node, simpacket);
    
    sim.wait();
    
    return sim;
}

/** Run the sub-system simulation consisting of 'nmoves' moves from 'moves'
    on the sub-system 'system', recording statistics if 'record_stats'
    is true. Run the simulation in the current thread */
SupraSubSim SupraSubSim::run(const SupraSubSystem &system, const SupraSubMoves &moves,
                             int nmoves, bool record_stats)
{
    return SupraSubSim::run( SupraSubSimPacket(system, moves, nmoves, record_stats) );
}
                    
/** Run the sub-system simulation consisting of 'nmoves' moves of 'move'
    on the sub-system 'system', recording statistics if 'record_stats'
    is true. Run the simulation in the current thread */
SupraSubSim SupraSubSim::run(const SupraSubSystem &system, const SupraSubMove &move,
                             int nmoves, bool record_stats)
{
    return SupraSubSim::run( SupraSubSimPacket(system, SameSupraSubMoves(move),
                                               nmoves, record_stats) );
}

/** Run the sub-system simulation consisting of 'nmoves' moves from 'moves'
    on the sub-system 'system', recording statistics if 'record_stats'
    is true. Run the simulation on the node 'node', returning a handle to
    the running simulation */
SupraSubSim SupraSubSim::run(Node &node,
                             const SupraSubSystem &system, const SupraSubMoves &moves,
                             int nmoves, bool record_stats)
{
    return SupraSubSim::run(node, SupraSubSimPacket(system,moves,nmoves,record_stats) );
}
                    
/** Run the sub-system simulation consisting of 'nmoves' moves of the move 'move'
    on the sub-system 'system', recording statistics if 'record_stats'
    is true. Run the simulation on the node 'node', returning a handle to
    the running simulation */
SupraSubSim SupraSubSim::run(Node &node,
                             const SupraSubSystem &system, const SupraSubMove &move,
                             int nmoves, bool record_stats)
{
    return SupraSubSim::run(node, SupraSubSimPacket(system, SameSupraSubMoves(move),
                                                    nmoves, record_stats) );
}

/** Abort the running simulation */
void SupraSubSim::abort()
{
    sim_promise.abort();
}

/** Stop the running simulation */
void SupraSubSim::stop()
{
    sim_promise.stop();
}

/** Wait for the simulation to complete */
void SupraSubSim::wait()
{
    sim_promise.wait();
}

/** Wait for the simulation to stop running, or for 'timeout'
    milliseconds to pass, whichever comes soonest. This returns
    whether or not the simulation has stopped */
bool SupraSubSim::wait(int timeout)
{
    return sim_promise.wait(timeout);
}

/** Return whether or not this simulation is running */
bool SupraSubSim::isRunning()
{
    return sim_promise.isRunning();
}

/** Return whether or not this simulation is in an error state */
bool SupraSubSim::isError()
{
    return sim_promise.isError();
}

/** Throw any error associated with this simulation - this does
    nothing if we are not in an error state */
void SupraSubSim::throwError()
{
    sim_promise.throwError();
}

/** Return whether or not the simulation was stopped */
bool SupraSubSim::wasStopped()
{
    return sim_promise.wasStopped();
}    

/** Return whether or not the simulation was aborted */
bool SupraSubSim::wasAborted()
{
    return sim_promise.wasAborted();
}

/** Return whether or not the simulation has finished
    (completed all of the moves) */
bool SupraSubSim::hasFinished()
{
    if (this->isRunning())
        return false;
        
    else
    {
        //we aren't running any more - lets see what happened
        if (this->isError() or this->wasAborted())
            return false;
            
        try
        {
            SupraSubSimPacket sim = this->result();
            
            return sim.nSubCompleted() == sim.nSubMoves();
        }
        catch(...)
        {
            return false;
        }
    }
}

/** Return the progress of the simulation (as a percentage) */
float SupraSubSim::progress()
{
    return sim_promise.progress();
}

/** Return the initial input simulation WorkPacket */
SupraSubSimPacket SupraSubSim::input()
{
    if (sim_promise.isNull())
        return SupraSubSimPacket();
        
    else
    {
        WorkPacket initial_packet = sim_promise.input();
        
        if (initial_packet.isNull())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the input simulation WorkPacket? How has "
                "it become null?"), CODELOC );
        }
        
        if (not initial_packet.isA<SupraSubSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the input simulation WorkPacket? How has "
                "it turned into a %1?").arg(initial_packet.base().what()),
                    CODELOC );
        }
    
        return initial_packet.asA<SupraSubSimPacket>();
    }
}

/** Return the simulation WorkPacket from an intermediate point along
    the simulation. This will throw an error if the simulation is in an
    error state, and the initial packet if the simulation 
    was aborted */
SupraSubSimPacket SupraSubSim::interimResult()
{
    if (sim_promise.isNull())
        return SupraSubSimPacket();
        
    else
    {
        WorkPacket interim_packet = sim_promise.interimResult();
        
        if (interim_packet.wasAborted())
        {
            return this->input();
        }
        else if (interim_packet.isError())
        {
            interim_packet.throwError();
        }
        
        if (interim_packet.isNull())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the interim simulation WorkPacket? How has "
                "it become null?"), CODELOC );
        }
        
        if (not interim_packet.isA<SupraSubSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the interim simulation WorkPacket? How has "
                "it turned into a %1?").arg(interim_packet.base().what()),
                    CODELOC );
        }
    
        return interim_packet.asA<SupraSubSimPacket>();
    }
}

/** Return the final result of the simulation. This blocks until
    the simulation has stopped, and will throw an exception if the
    simulation is in an error state. This returns the initial
    simulation WorkPacket if the simulation was aborted */
SupraSubSimPacket SupraSubSim::result()
{
    if (sim_promise.isNull())
        return SupraSubSimPacket();
        
    else
    {
        WorkPacket result_packet = sim_promise.result();
        
        if (result_packet.wasAborted())
        {
            return this->input();
        }
        else if (result_packet.isError())
        {
            result_packet.throwError();
        }
        
        if (result_packet.isNull())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the simulation result WorkPacket? How has "
                "it become null?"), CODELOC );
        }
        
        if (not result_packet.isA<SupraSubSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the simulation result WorkPacket? How has "
                "it turned into a %1?").arg(result_packet.base().what()),
                    CODELOC );
        }
    
        return result_packet.asA<SupraSubSimPacket>();
    }
}

/** Return the sub-system in the state it was in before the simulation started */
SupraSubSystemPtr SupraSubSim::initialSystem()
{
    return this->input().subSystem();
}

/** Return the Moves in the state they were in before the simulation started */
SupraSubMovesPtr SupraSubSim::initialMoves()
{
    return this->input().subMoves();
}

/** Return the current state of the System (updated while the simulation
    is running). This will throw an exception if the system hits an 
    error state */
SupraSubSystemPtr SupraSubSim::interimSystem()
{
    return this->interimResult().subSystem();
}

/** Return the current state of the moves (updated while the simulation
    is running). This will throw an exception if the system hits an 
    error state */
SupraSubMovesPtr SupraSubSim::interimMoves()
{
    return this->interimResult().subMoves();
}

/** Return the final state of the system after the simulation. This
    blocks until the simulation has finished and will throw an
    exception if the system hits an error state */
SupraSubSystemPtr SupraSubSim::system()
{
    return this->result().subSystem();
}

/** Return the final state of the moves after the simulation. This
    blocks until the simulation has finished and will throw an 
    exception if the system hits an error state */
SupraSubMovesPtr SupraSubSim::moves()
{
    return this->result().subMoves();
}
