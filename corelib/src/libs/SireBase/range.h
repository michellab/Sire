/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#ifndef SIREBASE_RANGE_H
#define SIREBASE_RANGE_H

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class Range;
}

QDataStream& operator<<(QDataStream&, const SireBase::Range&);
QDataStream& operator>>(QDataStream&, SireBase::Range&);

namespace SireBase
{

class SimpleRange;

typedef PropPtr<Range> RangePtr;

/** This class represents a range of integers.

    @author Christopher Woods
*/
class SIREBASE_EXPORT Range : public SireBase::Property
{
friend QDataStream& ::operator<<(QDataStream&, const SireBase::Range&);
friend QDataStream& ::operator>>(QDataStream&, SireBase::Range&);

public:
    Range();
    Range(const Range &other);

    virtual ~Range();
    
    virtual Range* clone() const=0;
    
    virtual qint64 next()=0;

    virtual bool atEnd() const=0;

    virtual RangePtr populate(int nitems) const=0;

    static SimpleRange null();

protected:
    Range& operator=(const Range &other);
};

}

SIRE_EXPOSE_CLASS( SireBase::Range )

SIRE_EXPOSE_PROPERTY( SireBase::RangePtr, SireBase::Range )

SIRE_END_HEADER

#endif
