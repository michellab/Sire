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

#include "ranges.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireID;
using namespace SireStream;

///////////////
/////////////// Implementation of SimpleRange
///////////////

static const RegisterMetaType<SimpleRange> r_simple;

QDataStream &operator<<(QDataStream &ds, const SimpleRange &range)
{
    writeHeader(ds, r_simple, 1);
    
    ds << range.strtval << range.endval << range.incr
       << static_cast<const Range&>(range);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, SimpleRange &range)
{
    VersionID v = readHeader(ds, r_simple);
    
    if (v == 1)
    {
        ds >> range.strtval >> range.endval >> range.incr
           >> static_cast<Range&>(range);
    }
    else
        throw version_error(v, "1", r_simple, CODELOC);
    
    return ds;
}

/** Construct a null range */
SimpleRange::SimpleRange()
            : ConcreteProperty<SimpleRange,Range>(),
              strtval(-1), endval(-1), incr(0)
{}

/** Construct a range that represents a single value, i */
SimpleRange::SimpleRange(qint64 i)
            : ConcreteProperty<SimpleRange,Range>(),
              strtval(i), endval(i+1), incr(1)
{}

/** Construct a range from [start,end) in units of increment */
SimpleRange::SimpleRange(qint64 start, qint64 end, qint64 increment)
            : ConcreteProperty<SimpleRange,Range>(),
              strtval(start), endval(end), incr(increment)
{
    if (increment == 0)
    {
        throw SireError::invalid_arg( QObject::tr(
                "You cannot construct a range using an increment of zero!"), CODELOC );
    }
}

/** Copy constructor */
SimpleRange::SimpleRange(const SimpleRange &other)
            : ConcreteProperty<SimpleRange,Range>(other),
              strtval(other.strtval), endval(other.endval), incr(other.incr)
{}

/** Destructor */
SimpleRange::~SimpleRange()
{}

/** Copy assignment operator */
SimpleRange& SimpleRange::operator=(const SimpleRange &other)
{
    if (this != &other)
    {
        strtval = other.strtval;
        endval = other.endval;
        incr = other.incr;
        Range::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool SimpleRange::operator==(const SimpleRange &other) const
{
    return strtval == other.strtval and
           endval == other.endval and
           incr == other.incr and
           Range::operator==(other);
}

/** Comparison operator */
bool SimpleRange::operator!=(const SimpleRange &other) const
{
    return not operator==(other);
}

/** Return a string representation */
QString SimpleRange::toString() const
{
    if (incr != 0)
    {
        if (strtval + incr == endval)
        {
            return QObject::tr("%1").arg(strtval);
        }
        else if (incr != 1)
        {
            return QObject::tr("[%1,%2,%3)").arg(strtval).arg(endval).arg(incr);
        }
        else
        {
            return QObject::tr("[%1,%2)").arg(strtval).arg(endval);
        }
    }
    else
        return QObject::tr("Range::null");
}

SimpleRange* SimpleRange::clone() const
{
    return new SimpleRange(*this);
}

const char* SimpleRange::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SimpleRange>() );
}

const char* SimpleRange::what() const
{
    return SimpleRange::typeName();
}

/** Return whether or not this range is at its end */
bool SimpleRange::atEnd() const
{
    return incr == 0;
}

/** Return the next value in the range */
qint64 SimpleRange::next()
{
    if (incr != 0)
    {
        qint64 val = strtval;

        strtval += incr;
        
        if (incr > 0)
        {
            if (strtval >= endval)
            {
                incr = 0;
            }
        }
        else
        {
            if (strtval <= endval)
            {
                incr = 0;
            }
        }
        
        return val;
    }
    else
        throw SireError::invalid_index( QObject::tr(
                "No more iterations available for this range!"), CODELOC );
    
    return 0;
}

/** Return a copy of this range that is populate with passed number of values */
RangePtr SimpleRange::populate(int nvalues) const
{
    if (incr == 0)
    {
        throw SireError::invalid_index( QObject::tr(
                "Cannot populate a null range!"), CODELOC );
    }

    SimpleRange ret(*this);
    ret.strtval = Index(strtval).map(nvalues);
    ret.endval = Index(endval).map(nvalues);

    if (incr > 0 and endval < 0)
        ret.endval += 1;

    return ret;
}
