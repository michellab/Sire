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

#include "range.h"
#include "ranges.h"

#include "SireStream/datastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Range> r_range( MAGIC_ONLY, "SireBase::Range" );

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const Range &range)
{
    writeHeader(ds, r_range, 1);
    
    ds << static_cast<const Property&>(range);
    
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, Range &range)
{
    VersionID v = readHeader(ds, r_range);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(range);
    }
    else
        throw version_error(v, "1", r_range, CODELOC);
    
    return ds;
}

/** Constructor */
Range::Range() : Property()
{}

/** Copy constructor */
Range::Range(const Range &other) : Property(other)
{}

/** Destructor */
Range::~Range()
{}

/** Copy assignment */
Range& Range::operator=(const Range &other)
{
    Property::operator=(other);
    return *this;
}

/** Return a null simple range for null */
const Range& Range::null()
{
    return *(create_shared_null<SimpleRange>());;
}

/** Return the range that represents the single value 'i' */
RangePtr Range::create(qint64 i)
{
    return SimpleRange(i);
}

/** Return the range that represents the range from [start,end) */
RangePtr Range::create(qint64 start, qint64 end)
{
    return SimpleRange(start,end);
}

/** Return the range that represents the range from [start,end,increment) */
RangePtr Range::create(qint64 start, qint64 end, qint64 increment)
{
    return SimpleRange(start,end,increment);
}
