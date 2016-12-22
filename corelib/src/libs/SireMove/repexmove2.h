/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#ifndef SIREMOVE_REPEXMOVE2_H
#define SIREMOVE_REPEXMOVE2_H

#include "SireMaths/rangenerator.h"

#include "SireMove/supramove.h"
#include "SireMove/suprasubmove.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class RepExMove2;
}

QDataStream& operator<<(QDataStream&, const SireMove::RepExMove2&);
QDataStream& operator>>(QDataStream&, SireMove::RepExMove2&);

namespace SireMove
{

class Replicas;
class Replica;

using SireMaths::RanGenerator;

/** This class is used to perform replica exchange moves on a collection
    of Replicas. Each move involves running a block of sampling
    on each of the replicas, and then performing replice exchange swaps
    and tests between pairs.
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT RepExMove2
        : public SireBase::ConcreteProperty<RepExMove2,SupraMove>
{

friend QDataStream& ::operator<<(QDataStream&, const RepExMove2&);
friend QDataStream& ::operator>>(QDataStream&, RepExMove2&);

public:
    RepExMove2();
    
    RepExMove2(const RepExMove2 &other);
    
    ~RepExMove2();
    
    RepExMove2& operator=(const RepExMove2 &other);

    bool operator==(const RepExMove2 &other) const;
    bool operator!=(const RepExMove2 &other) const;

    static const char* typeName();
    
    int nAttempted() const;
    int nAccepted() const;
    int nRejected() const;
    
    double acceptanceRatio() const;
    
    void clearStatistics();
    
    void setSwapMonitors(bool swap_monitors);
    
    bool swapMovesDisabled() const;
    void setDisableSwaps(bool disable);
    
    QString toString() const;
    
    void setGenerator(const RanGenerator &generator);
    const RanGenerator& generator() const;

    void move(SupraSystem &system, int nmoves, bool record_stats);

private:
    void performMove(Replicas &replicas, bool record_stats);

    /** The random number generator used to accept or reject the moves */
    RanGenerator rangenerator;
    
    /** The number of times a replica exchange move has been accepted */
    quint32 naccept;
    
    /** The number of times a replica exchange move has been rejected */
    quint32 nreject;

    /** Whether or not to swap the system monitors when we swap replicas
         - by default we leave the monitors with the systems */
    bool swap_monitors;
    
    /** Whether or not to disable RETI tests. This is useful when you want
        to just use this to RUN TI on a lot of replicas in parallel */
    bool disable_swaps;
};

} // end of namespace SireMove

Q_DECLARE_METATYPE( SireMove::RepExMove2 )

SIRE_EXPOSE_CLASS( SireMove::RepExMove2 )

#endif
