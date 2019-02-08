/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREMOVE_TITRATIONMOVE_H
#define SIREMOVE_TITRATIONMOVE_H

#include "montecarlo.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class TitrationMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::TitrationMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::TitrationMove&);

namespace SireMove
{

/** This class performs a Monte Carlo titration move. This moves
    a charge from one place to another by swapping the coordinates
    of once molecule with another, e.g. swapping a charge with a water.
    This allows ions to move quickly through a simulation box, and for
    ions to equilibrate between boxes (e.g. during a WSRC calcualtion)
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT TitrationMove
        : public SireBase::ConcreteProperty<TitrationMove,MonteCarlo>
{

friend SIREMOVE_EXPORT QDataStream& ::operator<<(QDataStream&, const TitrationMove&);
friend SIREMOVE_EXPORT QDataStream& ::operator>>(QDataStream&, TitrationMove&);

public:
    TitrationMove();
    TitrationMove(const TitrationMove &other);
    
    ~TitrationMove();

    static const char* typeName();
    
    const char* what() const;
    
    TitrationMove& operator=(const TitrationMove &other);
    
    bool operator==(const TitrationMove &other) const;
    bool operator!=(const TitrationMove &other) const;
    
    QString toString() const;
    
    void move(System &system, int nmoves, bool record_stats=true);

protected:
    void _pvt_setTemperature(const SireUnits::Dimension::Temperature &temperature);
};

} // end of namespace SireMove

Q_DECLARE_METATYPE( SireMove::TitrationMove )

SIRE_EXPOSE_CLASS( SireMove::TitrationMove )

SIRE_END_HEADER

#endif
