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

#ifndef SIREMOVE_SUPRAMOVES_H
#define SIREMOVE_SUPRAMOVES_H

#include "supramove.h"

SIRE_BEGIN_HEADER

namespace SireMove
{
class SupraMoves;
class SameSupraMoves;
}

QDataStream& operator<<(QDataStream&, const SireMove::SupraMoves&);
QDataStream& operator>>(QDataStream&, SireMove::SupraMoves&);

QDataStream& operator<<(QDataStream&, const SireMove::SameSupraMoves&);
QDataStream& operator>>(QDataStream&, SireMove::SameSupraMoves&);

namespace SireMove
{

/** This is the base class of all sets of moves that can be applied
    to supra-systems
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SupraMoves : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const SupraMoves&);
friend QDataStream& ::operator>>(QDataStream&, SupraMoves&);

public:
    SupraMoves();
    
    SupraMoves(const SupraMoves &other);
    
    virtual ~SupraMoves();
    
    virtual SupraMoves* clone() const=0;
    
    static const char* typeName()
    {
        return "SireMove::SupraMoves";
    }
    
    const SupraMove& operator[](int i) const;
    
    int count() const;
    int size() const;
    int nSubMoveTypes() const;
    
    virtual int nMoves() const=0;
    
    virtual QString toString() const=0;
    
    virtual void move(SupraSystem &system, int nmoves, bool record_stats=true)=0;
                      
    virtual void clearStatistics()=0;
    
    virtual QList<SupraMovePtr> moves() const=0;
    
    static const SameSupraMoves& null();

protected:
    SupraMoves& operator=(const SupraMoves &other);
    
    bool operator==(const SupraMoves &other) const;
    bool operator!=(const SupraMoves &other) const;
};

/** This SupraMoves object is used to apply the same SupraMove
    repeatedly to a SupraSystem
    
    @author Christopher Woods
*/
class SIREMOVE_EXPORT SameSupraMoves
            : public SireBase::ConcreteProperty<SameSupraMoves,SupraMoves>
{

friend QDataStream& ::operator<<(QDataStream&, const SameSupraMoves&);
friend QDataStream& ::operator>>(QDataStream&, SameSupraMoves&);

public:
    SameSupraMoves();
    SameSupraMoves(const SupraMove &move);
    
    SameSupraMoves(const SameSupraMoves &other);
    
    ~SameSupraMoves();
    
    SameSupraMoves& operator=(const SameSupraMoves &other);
    
    bool operator==(const SameSupraMoves &other) const;
    bool operator!=(const SameSupraMoves &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    void move(SupraSystem &system, int nmoves, bool record_stats=true);

    int nMoves() const;

    void clearStatistics();
    
    QList<SupraMovePtr> moves() const;

private:
    /** The move that will be repeatedly applied to the SupraSystem */
    SupraMovePtr mv;
};

typedef SireBase::PropPtr<SupraMoves> SupraMovesPtr;

}

Q_DECLARE_METATYPE( SireMove::SameSupraMoves )

SIRE_EXPOSE_CLASS( SireMove::SupraMoves )
SIRE_EXPOSE_CLASS( SireMove::SameSupraMoves )

SIRE_EXPOSE_PROPERTY( SireMove::SupraMovesPtr, SireMove::SupraMoves )

SIRE_END_HEADER

#endif
