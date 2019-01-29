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

#ifndef SIREMOVE_SUPRASUBMOVE_H
#define SIREMOVE_SUPRASUBMOVE_H

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraSubMove;
class NullSupraSubMove;
}

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::SupraSubMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::SupraSubMove&);

SIREMOVE_EXPORT QDataStream& operator<<(QDataStream&, const SireMove::NullSupraSubMove&);
SIREMOVE_EXPORT QDataStream& operator>>(QDataStream&, SireMove::NullSupraSubMove&);

namespace SireMove
{

class SupraSubSystem;

/** This is the base class of the controller for the sub-moves
    that are performed on the sub-systems of each SupraSystem
    as part of a SupraMove
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraSubMove : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const SupraSubMove&);
friend QDataStream& ::operator>>(QDataStream&, SupraSubMove&);

public:
    SupraSubMove();
    
    SupraSubMove(const SupraSubMove &other);
    
    virtual ~SupraSubMove();
    
    virtual SupraSubMove* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::SupraSubMove";
    }
    
    virtual QString toString() const=0;
    
    virtual void clearStatistics();
    
    virtual void move(SupraSubSystem &system, int n_supra_moves,
                      int n_supra_moves_per_block,
                      bool record_stats=true)=0;

    static const NullSupraSubMove& null();
    
protected:
    SupraSubMove& operator=(const SupraSubMove &other);
    
    bool operator==(const SupraSubMove &other) const;
    bool operator!=(const SupraSubMove &other) const; 

};

/** This is a null move that doesn't move a SupraSubSystem...

    @author Christopher Woods 
*/
class SIREMOVE_EXPORT NullSupraSubMove
        : public SireBase::ConcreteProperty<NullSupraSubMove,SupraSubMove>
{

friend QDataStream& ::operator<<(QDataStream&, const NullSupraSubMove&);
friend QDataStream& ::operator>>(QDataStream&, NullSupraSubMove&);

public:
    NullSupraSubMove();
    
    NullSupraSubMove(const NullSupraSubMove &other);
    
    ~NullSupraSubMove();
    
    NullSupraSubMove& operator=(const NullSupraSubMove &other);
    
    bool operator==(const NullSupraSubMove &other) const;
    bool operator!=(const NullSupraSubMove &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    void move(SupraSubSystem &system, int n_supra_moves, 
              int n_supra_moves_per_block, bool record_stats);
};

typedef SireBase::PropPtr<SupraSubMove> SupraSubMovePtr;

}

Q_DECLARE_METATYPE( SireMove::NullSupraSubMove )

SIRE_EXPOSE_CLASS( SireMove::SupraSubMove )
SIRE_EXPOSE_CLASS( SireMove::NullSupraSubMove )

SIRE_EXPOSE_PROPERTY( SireMove::SupraSubMovePtr, SireMove::SupraSubMove )

SIRE_END_HEADER

#endif
