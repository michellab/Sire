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

#ifndef SIREMOVE_SUPRAMOVE_H
#define SIREMOVE_SUPRAMOVE_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraMove;
class NullSupraMove;
}

QDataStream& operator<<(QDataStream&, const SireMove::SupraMove&);
QDataStream& operator>>(QDataStream&, SireMove::SupraMove&);

QDataStream& operator<<(QDataStream&, const SireMove::NullSupraMove&);
QDataStream& operator>>(QDataStream&, SireMove::NullSupraMove&);

namespace SireMove
{

class SupraSystem;

/** This is the base class of all supra-system moves (supra-moves).
    A supra-move is a move that is applied to a supra-system
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraMove : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const SupraMove&);
friend QDataStream& ::operator>>(QDataStream&, SupraMove&);

public:
    SupraMove();
    
    SupraMove(const SupraMove &other);
    
    virtual ~SupraMove();
    
    virtual SupraMove* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::SupraMove";
    }

    int nMoves() const;

    virtual QString toString() const=0;

    virtual void move(SupraSystem &system, int nmoves, 
                      bool record_stats=true)=0;

    virtual void clearStatistics();

    static const NullSupraMove& null();

protected:
    SupraMove& operator=(const SupraMove &other);
    
    bool operator==(const SupraMove &other) const;
    bool operator!=(const SupraMove &other) const;

    void incrementNMoves(int nmoves);

private:
    /** The total number of moves performed using this object */
    quint32 nmoves;
};

/** This is a null supra move, which does nothing

    @author Christopher Woods
*/
class SIREMOVE_EXPORT NullSupraMove
         : public SireBase::ConcreteProperty<NullSupraMove,SupraMove>
{

friend QDataStream& ::operator<<(QDataStream&, const NullSupraMove&);
friend QDataStream& ::operator>>(QDataStream&, NullSupraMove&);

public:
    NullSupraMove();
    
    NullSupraMove(const NullSupraMove &other);
    
    ~NullSupraMove();
    
    NullSupraMove& operator=(const NullSupraMove &other);
    
    bool operator==(const NullSupraMove &other) const;
    bool operator!=(const NullSupraMove &other) const;
    
    static const char* typeName();

    QString toString() const;

    void move(SupraSystem &system, int nmoves, bool record_stats=true);
};

typedef SireBase::PropPtr<SupraMove> SupraMovePtr;

}

Q_DECLARE_METATYPE( SireMove::NullSupraMove )

SIRE_EXPOSE_CLASS( SireMove::SupraMove )
SIRE_EXPOSE_CLASS( SireMove::NullSupraMove )

SIRE_EXPOSE_PROPERTY( SireMove::SupraMovePtr, SireMove::SupraMove )

SIRE_END_HEADER

#endif

