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

#ifndef SIREBASE_RANGES_H
#define SIREBASE_RANGES_H

#include "range.h"

namespace SireBase
{
class SimpleRange;
class SetRange;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::SimpleRange&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::SimpleRange&);

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::SetRange&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::SetRange&);

namespace SireBase
{

/** This class represents a simple range from start to end in steps
    of increment
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT SimpleRange : public SireBase::ConcreteProperty<SimpleRange,Range>
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const SireBase::SimpleRange&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, SireBase::SimpleRange&);

public:
    SimpleRange();
    SimpleRange(qint64 i);
    SimpleRange(qint64 start, qint64 end, qint64 increment=1);
    
    SimpleRange(const SimpleRange &other);
    
    ~SimpleRange();
    
    SimpleRange& operator=(const SimpleRange &other);
    
    bool operator==(const SimpleRange &other) const;
    bool operator!=(const SimpleRange &other) const;
    
    QString toString() const;
    
    SimpleRange* clone() const;
    
    static const char* typeName();
    const char* what() const;
    
    qint64 next();
    
    bool atEnd() const;
    
    RangePtr populate(int nvalues) const;
    
private:
    qint64 strtval, endval, incr;
};

}

Q_DECLARE_METATYPE( SireBase::SimpleRange )

SIRE_EXPOSE_CLASS( SireBase::SimpleRange )

#endif
