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

#ifndef SIREMOVE_SUPRASIMPACKET_H
#define SIREMOVE_SUPRASIMPACKET_H

#include "SireCluster/workpacket.h"

#include "suprasystem.h"
#include "supramoves.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraSimPacket;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::SupraSimPacket&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::SupraSimPacket&);

namespace SireMove
{

/** This is a workpacket that is used to run part of a supra-simulation

    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSimPacket : public SireCluster::WorkPacketBase
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const SupraSimPacket&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, SupraSimPacket&);

public:
    SupraSimPacket();
    
    SupraSimPacket(const SupraSystem &suprasystem,
                   const SupraMoves &supramoves,
                   int nmoves, bool record_stats);
                   
    SupraSimPacket(const SupraSimPacket &other);
    
    ~SupraSimPacket();
    
    SupraSimPacket& operator=(const SupraSimPacket &other);
    
    bool operator==(const SupraSimPacket &other) const;
    bool operator!=(const SupraSimPacket &other) const;
    
    static const char* typeName();
    
    const char* what() const
    {
        return SupraSimPacket::typeName();
    }
    
    SupraSimPacket* clone() const;
    
    bool shouldPack() const;
    int approximatePacketSize() const;

    const SupraSystem& system() const;
    const SupraMoves& moves() const;
    
    int nMoves() const;
    int nCompleted() const;

    bool recordingStatistics() const;
    
    bool hasFinished() const;

protected:
    float chunk();
    
private:
    /** The supra-system being simulated */
    SupraSystemPtr supra_system;
    
    /** The moves applied to the supra-system */
    SupraMovesPtr supra_moves;
    
    /** The number of moves to be run on the supra-system */
    quint32 n_supra_moves;
    
    /** The number of supra-moves already run on the supra-system */
    quint32 ncompleted;
    
    /** Whether or not to record sub statistics */
    bool record_stats;
};

}

Q_DECLARE_METATYPE( SireMove::SupraSimPacket )

SIRE_EXPOSE_CLASS( SireMove::SupraSimPacket )

SIRE_END_HEADER

#endif
