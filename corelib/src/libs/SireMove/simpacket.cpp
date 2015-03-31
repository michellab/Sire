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

#include "simpacket.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireSystem;
using namespace SireCluster;
using namespace SireStream;

static const RegisterMetaType<SimPacket> r_simpacket;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const SimPacket &simpacket)
{
    writeHeader(ds, r_simpacket, 2);
    
    SharedDataStream sds(ds);
    
    if (simpacket.sim_store_was_packed)
    {
        SimStore packed_store = simpacket.sim_store;
        packed_store.pack();
        
        sds << packed_store;
    }
    else
        sds << simpacket.sim_store;
    
    sds << simpacket.nmoves << simpacket.ncompleted
        << simpacket.nmoves_per_chunk << simpacket.record_stats;
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, SimPacket &simpacket)
{
    VersionID v = readHeader(ds, r_simpacket);
    
    if (v == 2)
    {
        SharedDataStream sds(ds);
        
        sds >> simpacket.sim_store
            >> simpacket.nmoves >> simpacket.ncompleted
            >> simpacket.nmoves_per_chunk >> simpacket.record_stats;
    }
    else if (v == 1)
    {
        SharedDataStream sds(ds);

        System sim_system;
        MovesPtr sim_moves;
        
        sds >> sim_system >> sim_moves
            >> simpacket.nmoves >> simpacket.ncompleted
            >> simpacket.nmoves_per_chunk >> simpacket.record_stats;
            
        simpacket.sim_store = SimStore(sim_system, sim_moves);
    }
    else
        throw version_error( v, "1", r_simpacket, CODELOC );
        
    return ds;
}

/** Null constructor */
SimPacket::SimPacket() 
          : WorkPacketBase(), nmoves(0), ncompleted(0),
            nmoves_per_chunk(0), record_stats(true),
            sim_store_was_packed(false)
{}

/** Construct a workpacket that runs 'nmoves' of the Moves 'moves' on the 
    passed System 'system', optionally recording simulation statistics
    if 'record_stats' is true */
SimPacket::SimPacket(const System &system, const Moves &moves,
                     int n_moves, bool recording_stats)
          : WorkPacketBase(), sim_store(system,moves),
            ncompleted(0), nmoves_per_chunk(100), record_stats(recording_stats),
            sim_store_was_packed(false)
{
    if (n_moves > 0)
        nmoves = n_moves;
    else
        nmoves = 0;
}
          
/** Construct a workpacket that runs 'nmoves' of the Moves 'moves' on the 
    passed System 'system', optionally recording simulation statistics
    if 'record_stats' is true, and running 'nmoves_per_chunk' moves
    for each chunk */
SimPacket::SimPacket(const System &system, const Moves &moves,
                     int n_moves, int n_moves_per_chunk, bool recording_stats)
          : WorkPacketBase(), sim_store(system,moves),
            ncompleted(0), record_stats(recording_stats),
            sim_store_was_packed(false)
{
    if (n_moves > 0)
        nmoves = n_moves;
    else
        nmoves = 0;

    if (n_moves_per_chunk > 0)
        nmoves_per_chunk = n_moves_per_chunk;
    else
        nmoves_per_chunk = 1;
}

/** Construct a workpacket that runs 'nmoves' of the Moves on the 
    System, both contained in 'simstore', optionally recording simulation statistics
    if 'record_stats' is true */
SimPacket::SimPacket(const SimStore &simstore,
                     int n_moves, bool recording_stats)
          : WorkPacketBase(), sim_store(simstore),
            ncompleted(0), nmoves_per_chunk(100), record_stats(recording_stats),
            sim_store_was_packed(false)
{
    if (n_moves > 0)
        nmoves = n_moves;
    else
        nmoves = 0;
}
          
/** Construct a workpacket that runs 'nmoves' of the Moves  on the 
    System, both contained in 'simstore', optionally recording simulation statistics
    if 'record_stats' is true, and running 'nmoves_per_chunk' moves
    for each chunk */
SimPacket::SimPacket(const SimStore &simstore,
                     int n_moves, int n_moves_per_chunk, bool recording_stats)
          : WorkPacketBase(), sim_store(simstore),
            ncompleted(0), record_stats(recording_stats),
            sim_store_was_packed(false)
{
    if (n_moves > 0)
        nmoves = n_moves;
    else
        nmoves = 0;

    if (n_moves_per_chunk > 0)
        nmoves_per_chunk = n_moves_per_chunk;
    else
        nmoves_per_chunk = 1;
}

/** Copy constructor */
SimPacket::SimPacket(const SimPacket &other)
          : WorkPacketBase(other), sim_store(other.sim_store),
            nmoves(other.nmoves),
            ncompleted(other.ncompleted), nmoves_per_chunk(other.nmoves_per_chunk),
            record_stats(other.record_stats),
            sim_store_was_packed(other.sim_store_was_packed)
{}

/** Destructor */
SimPacket::~SimPacket()
{}

/** Copy assignment operator */
SimPacket& SimPacket::operator=(const SimPacket &other)
{
    if (this != &other)
    {
        sim_store = other.sim_store;
        nmoves = other.nmoves;
        ncompleted = other.ncompleted;
        nmoves_per_chunk = other.nmoves_per_chunk;
        record_stats = other.record_stats;
        sim_store_was_packed = other.sim_store_was_packed;
        
        WorkPacketBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SimPacket::operator==(const SimPacket &other) const
{
    return this == &other or
           (sim_store == other.sim_store and
            nmoves == other.nmoves and ncompleted == other.ncompleted and
            nmoves_per_chunk == other.nmoves_per_chunk and
            record_stats == other.record_stats and
            sim_store_was_packed == other.sim_store_was_packed);
}

/** Comparison operator */
bool SimPacket::operator!=(const SimPacket &other) const
{
    return not this->operator==(other);
}

/** Only compress this workpacket if the SimStore is not already packed */
bool SimPacket::shouldPack() const
{
    return not (sim_store_was_packed or sim_store.isPacked());
}

/** Because it takes too long to calculate the size of this
    packet, we say that it will be 32 MB - this is enough for 
    most cases */
int SimPacket::approximatePacketSize() const
{
    return 32 * 1024 * 1024;
}

/** Return the system being simulated */
System SimPacket::system() const
{
    return sim_store.system();
}

/** Return the moves being applied to the system */
MovesPtr SimPacket::moves() const
{
    return sim_store.moves();
}

/** Return both the system and moves together */
SimStore SimPacket::systemAndMoves() const
{
    return sim_store;
}

/** Return the number of moves being applied to the system */
int SimPacket::nMoves() const
{
    return nmoves;
}

/** Return the number of moves already run on the system */
int SimPacket::nCompleted() const
{
    return ncompleted;
}

/** Return the number of moves to apply for each chunk */
int SimPacket::nMovesPerChunk() const
{
    return nmoves_per_chunk;
}

/** Return whether or not simulation statistics will be recorded
    during the moves */
bool SimPacket::recordingStatistics() const
{
    return record_stats;
}

/** Return whether or not this simulation has finished */
bool SimPacket::hasFinished() const
{
    return ncompleted >= nmoves;
}

/** Run a chunk of the simulation, and return the progress */
float SimPacket::chunk()
{
    if (nmoves == 0)
        return 100;

    int n_to_run = qMin( nmoves-ncompleted, nmoves_per_chunk );

    if (n_to_run > 0)
    {
        if (sim_store.isPacked())
        {
            sim_store_was_packed = true;
        
            //extract the system and moves from the store
            sim_store.unpack();
        }
        
        System sim_system = sim_store.system();
        MovesPtr sim_moves = sim_store.moves();
    
        //run a chunk of moves
        sim_system = sim_moves.edit().move(sim_system, n_to_run, record_stats);
        
        sim_store.setSystemAndMoves(sim_system, sim_moves);
        
        //it all completed successfully :-)
        ncompleted += n_to_run;

        //we leave the store unpacked, as there is no point repacking
        //it between chunks
    }

    if (ncompleted >= nmoves)
    {
        //we have finished all of the moves, so repack the simstore
        //if necessary
        if (sim_store_was_packed)
        {
            sim_store.pack();
            sim_store_was_packed = false;
        }
    }

    return 100.0 * ( float(ncompleted) / float(nmoves) );
}

const char* SimPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SimPacket>() );
}

SimPacket* SimPacket::clone() const
{
    return new SimPacket(*this);
}
