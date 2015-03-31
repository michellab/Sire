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

#include "suprasubsimpacket.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMove;
using namespace SireCluster;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<SupraSubSimPacket> r_suprasubsimpacket;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const SupraSubSimPacket &suprasubsimpacket)
{
    writeHeader(ds, r_suprasubsimpacket, 1);
    
    SharedDataStream sds(ds);
    
    if (suprasubsimpacket.sub_system_was_packed)
    {
        SupraSubSystemPtr packed_sub_system = suprasubsimpacket.sub_system;
        packed_sub_system.edit().pack();

        sds << packed_sub_system;
    }
    else
        sds << suprasubsimpacket.sub_system;
    
    sds << suprasubsimpacket.sub_moves
        << suprasubsimpacket.n_sub_moves << suprasubsimpacket.ncompleted
        << suprasubsimpacket.record_stats
        << static_cast<const WorkPacketBase&>(suprasubsimpacket);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        SupraSubSimPacket &suprasubsimpacket)
{
    VersionID v = readHeader(ds, r_suprasubsimpacket);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> suprasubsimpacket.sub_system >> suprasubsimpacket.sub_moves
            >> suprasubsimpacket.n_sub_moves >> suprasubsimpacket.ncompleted
            >> suprasubsimpacket.record_stats
            >> static_cast<WorkPacketBase&>(suprasubsimpacket);
            
        suprasubsimpacket.sub_system_was_packed = false;
    }
    else
        throw version_error(v, "1", r_suprasubsimpacket, CODELOC);
        
    return ds;
}

/** Constructor */
SupraSubSimPacket::SupraSubSimPacket() 
                  : WorkPacketBase(), n_sub_moves(0), ncompleted(0), 
                    record_stats(false), sub_system_was_packed(false)
{}

/** Construct a work packet to perform 'nmoves' sub-moves (in 'moves') on 
    the sub-system 'system', recording statistics if 'record_stats' is true */
SupraSubSimPacket::SupraSubSimPacket(const SupraSubSystem &system,
                                     const SupraSubMoves &moves,
                                     int nmoves, bool record_statistics)
                  : WorkPacketBase(),
                    sub_system(system), sub_moves(moves),
                    n_sub_moves(nmoves), ncompleted(0), 
                    record_stats(record_statistics), sub_system_was_packed(false)
{}
  
/** Copy constructor */                
SupraSubSimPacket::SupraSubSimPacket(const SupraSubSimPacket &other)
                  : WorkPacketBase(other),
                    sub_system(other.sub_system), sub_moves(other.sub_moves),
                    n_sub_moves(other.n_sub_moves), ncompleted(other.ncompleted),
                    record_stats(other.record_stats),
                    sub_system_was_packed(other.sub_system_was_packed)
{}

/** Destructor */
SupraSubSimPacket::~SupraSubSimPacket()
{}

/** Copy assignment operator */
SupraSubSimPacket& SupraSubSimPacket::operator=(const SupraSubSimPacket &other)
{
    if (this != &other)
    {
        sub_system = other.sub_system;
        sub_moves = other.sub_moves;
        n_sub_moves = other.n_sub_moves;
        ncompleted = other.ncompleted;
        record_stats = other.record_stats;
        sub_system_was_packed = other.sub_system_was_packed;
        
        WorkPacketBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SupraSubSimPacket::operator==(const SupraSubSimPacket &other) const
{
    return (this == &other) or
           (sub_system == other.sub_system and
            sub_moves == other.sub_moves and
            n_sub_moves == other.n_sub_moves and
            ncompleted == other.ncompleted and
            record_stats == other.record_stats);
}

/** Comparison operator */
bool SupraSubSimPacket::operator!=(const SupraSubSimPacket &other) const
{
    return not this->operator==(other);
}

/** We probably shouldn't pack this workpacket as it is already
    heavily packed, and if it is large, then it is likely already
    packed to disk */
bool SupraSubSimPacket::shouldPack() const
{
    return false;
}

/** This is likely to be large */
int SupraSubSimPacket::approximatePacketSize() const
{
    //16 MB
    return 16*1024*1024;
}

/** Return the sub-system */
const SupraSubSystem& SupraSubSimPacket::subSystem() const
{
    return sub_system;
}

/** Return the moves */
const SupraSubMoves& SupraSubSimPacket::subMoves() const
{
    return sub_moves;
}

/** Return the number of sub-moves to be applied to the sub-system */
int SupraSubSimPacket::nSubMoves() const
{
    return n_sub_moves;
}

/** Return the number of completed sub-moves */
int SupraSubSimPacket::nSubCompleted() const
{
    return ncompleted;
}

/** Return whether or not we are recording statistics during the sub-moves */
bool SupraSubSimPacket::recordingSubStatistics() const
{
    return record_stats;
}

/** Return whether or not this work packet has finished */
bool SupraSubSimPacket::hasFinished() const
{
    return ncompleted >= n_sub_moves;
}

/** Perform a chunk of simulation - this performs just one
    sub-move, as the sub-move will itself consist of several
    moves on the system */
float SupraSubSimPacket::chunk()
{
    if (ncompleted >= n_sub_moves)
        return 100.0;

    //see if we need to unpack the system before running the moves
    if (sub_system->isPacked())
    {
        sub_system_was_packed = true;
        sub_system.edit().unpack();
    }
        
    sub_moves.edit().move( sub_system.edit(), 1, n_sub_moves - ncompleted,
                           record_stats );
    
    ++ncompleted;
    
    if (ncompleted >= n_sub_moves)
    {
        if (sub_system_was_packed)
        {
            sub_system.edit().pack();
            sub_system_was_packed = false;
        }
    }
    
    return 100.0 * ( float(ncompleted) / float(n_sub_moves) );
}

const char* SupraSubSimPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SupraSubSimPacket>() );
}

SupraSubSimPacket* SupraSubSimPacket::clone() const
{
    return new SupraSubSimPacket(*this);
}
