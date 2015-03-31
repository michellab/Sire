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

#ifndef SIREMOVE_SUPRASUBSIMPACKET_H
#define SIREMOVE_SUPRASUBSIMPACKET_H

#include "SireCluster/workpacket.h"

#include "suprasubsystem.h"
#include "suprasubmoves.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraSubSimPacket;
}

QDataStream& operator<<(QDataStream&, const SireMove::SupraSubSimPacket&);
QDataStream& operator>>(QDataStream&, SireMove::SupraSubSimPacket&);

namespace SireMove
{

/** This is a workpacket that is used to run part of a SupraSubSim simulation

    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSubSimPacket : public SireCluster::WorkPacketBase
{

friend QDataStream& ::operator<<(QDataStream&, const SupraSubSimPacket&);
friend QDataStream& ::operator>>(QDataStream&, SupraSubSimPacket&);

public:
    SupraSubSimPacket();
    
    SupraSubSimPacket(const SupraSubSystem &system,
                      const SupraSubMoves &moves,
                      int nmoves, bool record_stats=true);
                      
    SupraSubSimPacket(const SupraSubSimPacket &other);
    
    ~SupraSubSimPacket();
    
    SupraSubSimPacket& operator=(const SupraSubSimPacket &other);
    
    bool operator==(const SupraSubSimPacket &other) const;
    bool operator!=(const SupraSubSimPacket &other) const;
    
    static const char* typeName();
    
    const char* what() const
    {
        return SupraSubSimPacket::typeName();
    }
    
    SupraSubSimPacket* clone() const;
    
    bool shouldPack() const;
    int approximatePacketSize() const;

    const SupraSubSystem& subSystem() const;
    const SupraSubMoves& subMoves() const;
    
    int nSubMoves() const;
    int nSubCompleted() const;

    bool recordingSubStatistics() const;
    
    bool hasFinished() const;

protected:
    float chunk();
    
private:
    /** The subsystem being simulated */
    SupraSubSystemPtr sub_system;
    
    /** The moves applied to the subsystem */
    SupraSubMovesPtr sub_moves;
    
    /** The number of submoves to be run on the system */
    quint32 n_sub_moves;
    
    /** The number of submoves already run on the system */
    quint32 ncompleted;
    
    /** Whether or not to record sub statistics */
    bool record_stats;

    /** Whether or not the sub-system is packed */
    bool sub_system_was_packed;
};

}

Q_DECLARE_METATYPE( SireMove::SupraSubSimPacket )

SIRE_EXPOSE_CLASS( SireMove::SupraSubSimPacket )

SIRE_END_HEADER

#endif
