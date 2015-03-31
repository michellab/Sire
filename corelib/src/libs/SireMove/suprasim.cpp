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

#include "suprasim.h"
#include "supramoves.h"

#include "SireCluster/nodes.h"

#include "SireError/errors.h"

using namespace SireMove;
using namespace SireCluster;

/** Null constructor */
SupraSim::SupraSim()
{}

/** Internal constructor used to connect to a running simulation */
SupraSim::SupraSim(const SireCluster::Promise &promise)
         : sim_promise(promise)
{}

/** Copy constructor */
SupraSim::SupraSim(const SupraSim &other)
         : sim_promise(other.sim_promise)
{}

/** Destructor */
SupraSim::~SupraSim()
{}

/** Copy assignment operator */
SupraSim& SupraSim::operator=(const SupraSim &other)
{
    sim_promise = other.sim_promise;
    return *this;
}

/** Comparison operator */
bool SupraSim::operator==(const SupraSim &other) const
{
    return sim_promise == other.sim_promise;
}

/** Comparison operator */
bool SupraSim::operator!=(const SupraSim &other) const
{
    return sim_promise != other.sim_promise;
}

/** Run the supra-system simulation described in 'simpacket' on the
    node 'node', returning a handle to the running simulation */
SupraSim SupraSim::run(Node &node, const SupraSimPacket &simpacket)
{
    return SupraSim( node.startJob(simpacket) );
}

/** Run the supra-system simulation described in 'simpacket' in the 
    current thread */
SupraSim SupraSim::run(const SupraSimPacket &simpacket)
{
    Nodes nodes;
    
    ThisThread this_thread = nodes.borrowThisThread();
    
    if (nodes.isEmpty())
        throw SireError::unavailable_resource( QObject::tr(
            "This thread is unavailable for running a simulation. It is already "
            "busy doing something else!"), CODELOC );
            
    Node node = nodes.getNode();
    
    SupraSim sim = SupraSim::run(node, simpacket);
    
    sim.wait();
    
    return sim;
}

/** Run the supra-system simulation applying 'nmoves' moves from 'moves'
    on the supra-system 'system', recording statistics if 'record_stats' 
    is true. The simulation is run in the current thread */
SupraSim SupraSim::run(const SupraSystem &system, const SupraMoves &moves,
                       int nmoves, bool record_stats)
{
    return SupraSim::run( SupraSimPacket(system,moves,nmoves,record_stats) );
}

/** Run the supra-system simulation consisting of 'nmoves' move of the move 'move' 
    on the supra-system 'system', recording statistics if 'record_stats' 
    is true. The simulation is run in the current thread */
SupraSim SupraSim::run(const SupraSystem &system, const SupraMove &move,
                       int nmoves, bool record_stats)
{
    return SupraSim::run( SupraSimPacket(system,SameSupraMoves(move),
                                         nmoves,record_stats) );
}
                
/** Run the supra-system simulation applying 'nmoves' moves from 'moves'
    on the supra-system 'system', recording statistics if 'record_stats' 
    is true. The simulation is run on the node 'node' and a handle is 
    returned to the running simulation */
SupraSim SupraSim::run(Node &node,
                       const SupraSystem &system, const SupraMoves &moves,
                       int nmoves, bool record_stats)
{
    return SupraSim::run(node, SupraSimPacket(system,moves,nmoves,record_stats) );
}
                
/** Run the supra-system simulation consisting of 'nmoves' move of the move 'move' 
    on the supra-system 'system', recording statistics if 'record_stats' 
    is true. The simulation is run on the node 'node' and a handle is
    returned to the running simulation */
SupraSim SupraSim::run(Node &node,
                       const SupraSystem &system, const SupraMove &move,
                       int nmoves, bool record_stats)
{
    return SupraSim::run(node, SupraSimPacket(system, SameSupraMoves(move),
                                              nmoves, record_stats) );
}

/** Abort the simulation */
void SupraSim::abort()
{
    sim_promise.abort();
}

/** Stop the simulation */
void SupraSim::stop()
{
    sim_promise.stop();
}

/** Wait until the simulation has finished */
void SupraSim::wait()
{
    sim_promise.wait();
}

/** Wait for the simulation to stop running, or for 'timeout'
    milliseconds to pass, whichever comes soonest. This returns
    whether or not the simulation has stopped */
bool SupraSim::wait(int timeout)
{
    return sim_promise.wait(timeout);
}

/** Return whether or not this simulation is running */
bool SupraSim::isRunning()
{
    return sim_promise.isRunning();
}

/** Return whether or not this simulation is in an error state */
bool SupraSim::isError()
{
    return sim_promise.isError();
}

/** Throw any error associated with this simulation - this does
    nothing if we are not in an error state */
void SupraSim::throwError()
{
    sim_promise.throwError();
}

/** Return whether or not the simulation was stopped */
bool SupraSim::wasStopped()
{
    return sim_promise.wasStopped();
}

/** Return whether or not the simulation was aborted */
bool SupraSim::wasAborted()
{
    return sim_promise.wasAborted();
}

/** Return whether or not the simulation has finished
    (completed all of the moves) */
bool SupraSim::hasFinished()
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
            SupraSimPacket sim = this->result();
            
            return sim.nCompleted() == sim.nMoves();
        }
        catch(...)
        {
            return false;
        }
    }
}

/** Return the progress of the simulation (as a percentage) */
float SupraSim::progress()
{
    return sim_promise.progress();
}

/** Return the initial input simulation WorkPacket */
SupraSimPacket SupraSim::input()
{
    if (sim_promise.isNull())
        return SupraSimPacket();
        
    else
    {
        WorkPacket initial_packet = sim_promise.input();
        
        if (initial_packet.isNull())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the input simulation WorkPacket? How has "
                "it become null?"), CODELOC );
        }
        
        if (not initial_packet.isA<SupraSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the input simulation WorkPacket? How has "
                "it turned into a %1?").arg(initial_packet.base().what()),
                    CODELOC );
        }
    
        return initial_packet.asA<SupraSimPacket>();
    }
}

/** Return the simulation WorkPacket from an intermediate point along
    the simulation. This will throw an error if the simulation is in an
    error state, and the initial packet if the simulation 
    was aborted */
SupraSimPacket SupraSim::interimResult()
{
    if (sim_promise.isNull())
        return SupraSimPacket();
        
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
        
        if (not interim_packet.isA<SupraSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the interim simulation WorkPacket? How has "
                "it turned into a %1?").arg(interim_packet.base().what()),
                    CODELOC );
        }
    
        return interim_packet.asA<SupraSimPacket>();
    }
}

/** Return the final result of the simulation. This blocks until
    the simulation has stopped, and will throw an exception if the
    simulation is in an error state. This returns the initial
    simulation WorkPacket if the simulation was aborted */
SupraSimPacket SupraSim::result()
{
    if (sim_promise.isNull())
        return SupraSimPacket();
        
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
        
        if (not result_packet.isA<SupraSimPacket>())
        {
            throw SireError::program_bug( QObject::tr(
                "How could we lose the simulation result WorkPacket? How has "
                "it turned into a %1?").arg(result_packet.base().what()),
                    CODELOC );
        }
    
        return result_packet.asA<SupraSimPacket>();
    }
}

/** Return the System in the state it was in before the simulation started */
SupraSystemPtr SupraSim::initialSystem()
{
    return this->input().system();
}

/** Return the Moves in the state they were in before the simulation started */
SupraMovesPtr SupraSim::initialMoves()
{
    return this->input().moves();
}

/** Return the current state of the System (updated while the simulation
    is running). This will throw an exception if the system hits an 
    error state */
SupraSystemPtr SupraSim::interimSystem()
{
    return this->interimResult().system();
}

/** Return the current state of the moves (updated while the simulation
    is running). This will throw an exception if the system hits an 
    error state */
SupraMovesPtr SupraSim::interimMoves()
{
    return this->interimResult().moves();
}

/** Return the final state of the system after the simulation. This
    blocks until the simulation has finished and will throw an
    exception if the system hits an error state */
SupraSystemPtr SupraSim::system()
{
    return this->result().system();
}

/** Return the final state of the moves after the simulation. This
    blocks until the simulation has finished and will throw an 
    exception if the system hits an error state */
SupraMovesPtr SupraSim::moves()
{
    return this->result().moves();
}
