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

#include "suprasimpacket.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMove;
using namespace SireCluster;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<SupraSimPacket> r_suprasimpacket;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, 
                                        const SupraSimPacket &suprasimpacket)
{
    writeHeader(ds, r_suprasimpacket, 1);
    
    SharedDataStream sds(ds);
    
    sds << suprasimpacket.supra_system << suprasimpacket.supra_moves
        << suprasimpacket.n_supra_moves << suprasimpacket.ncompleted
        << suprasimpacket.record_stats
        << static_cast<const WorkPacketBase&>(suprasimpacket);
        
    return ds; 
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, SupraSimPacket &suprasimpacket)
{
    VersionID v = readHeader(ds, r_suprasimpacket);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> suprasimpacket.supra_system >> suprasimpacket.supra_moves
            >> suprasimpacket.n_supra_moves >> suprasimpacket.ncompleted
            >> suprasimpacket.record_stats
            >> static_cast<WorkPacketBase&>(suprasimpacket);
    }
    else
        throw version_error(v, "1", r_suprasimpacket, CODELOC);
        
    return ds;
}

/** Null constructor */
SupraSimPacket::SupraSimPacket() 
               : WorkPacketBase(), n_supra_moves(0), ncompleted(0), record_stats(false)
{}

/** Construct to perform the passed supra-moves on the passed supra-system.
    The packet is to run 'nmoves' blocks of these moves, recording
    statistics if 'record_stats' is true */
SupraSimPacket::SupraSimPacket(const SupraSystem &suprasystem,
                               const SupraMoves &supramoves,
                               int nmoves, bool recording_stats)
               : WorkPacketBase(), 
                 supra_system(suprasystem), supra_moves(supramoves),
                 n_supra_moves(nmoves), ncompleted(0), record_stats(recording_stats)
{}

/** Copy constructor */               
SupraSimPacket::SupraSimPacket(const SupraSimPacket &other)
               : WorkPacketBase(other),
                 supra_system(other.supra_system), 
                 supra_moves(other.supra_moves),
                 n_supra_moves(other.n_supra_moves),
                 ncompleted(other.ncompleted),
                 record_stats(other.record_stats)
{}

/** Destructor */
SupraSimPacket::~SupraSimPacket()
{}

/** Copy assignment operator */
SupraSimPacket& SupraSimPacket::operator=(const SupraSimPacket &other)
{
    if (this != &other)
    {
        supra_system = other.supra_system;
        supra_moves = other.supra_moves;
        n_supra_moves = other.n_supra_moves;
        ncompleted = other.ncompleted;
        record_stats = other.record_stats;
        
        WorkPacketBase::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SupraSimPacket::operator==(const SupraSimPacket &other) const
{
    return (this == &other) or
           (supra_system == other.supra_system and
            supra_moves == other.supra_moves and
            n_supra_moves == other.n_supra_moves and
            ncompleted == other.ncompleted and
            record_stats == other.record_stats);
}

/** Comparison operator */
bool SupraSimPacket::operator!=(const SupraSimPacket &other) const
{
    return not this->operator==(other);
}

/** This probably shouldn't be packed to disk, as there will be a lot 
    of data sharing between this packet and other copies at different
    stages of the simulation (I think) - it is also already heavily
    packed (e.g. to disk) and I don't want that data to be pulled
    into memory */
bool SupraSimPacket::shouldPack() const
{
    return false;
}

/** This will be large...! */
int SupraSimPacket::approximatePacketSize() const
{
    //128 MB
    return 128*1024*1024;
}

/** Return the supra-system being simulated */
const SupraSystem& SupraSimPacket::system() const
{
    return supra_system;
}

/** Return the supra-moves being applied to the supra-system */
const SupraMoves& SupraSimPacket::moves() const
{
    return supra_moves;
}

/** Return the number of supra-moves to be applied to the supra-system */
int SupraSimPacket::nMoves() const
{
    return n_supra_moves;
}

/** Return the number of supra-moves that have been completed so far */
int SupraSimPacket::nCompleted() const
{
    return ncompleted;
}

/** Set whether or not statistics are being recorded during the moves */
bool SupraSimPacket::recordingStatistics() const
{
    return record_stats;
}

/** Return whether or not the simulation has finished */
bool SupraSimPacket::hasFinished() const
{
    return ncompleted >= n_supra_moves;
}

/** Perform a chunk of simulation - this performs just one supra-move
    (as each move will contain lots of sub-moves) */
float SupraSimPacket::chunk()
{
    if (ncompleted >= n_supra_moves)
        return 100.0;
        
    supra_moves.edit().move( supra_system.edit(), 1, record_stats );
    
    ++ncompleted;
    
    return 100.0 * ( float(ncompleted) / float(n_supra_moves) );
}

const char* SupraSimPacket::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SupraSimPacket>() );
}

SupraSimPacket* SupraSimPacket::clone() const
{
    return new SupraSimPacket(*this);
}
