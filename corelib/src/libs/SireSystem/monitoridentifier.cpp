/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "monitoridentifier.h"
#include "monitoridx.h"
#include "monitorname.h"

#include "systemmonitors.h"

#include "SireSystem/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/streampolypointer.hpp"

using namespace SireSystem;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<MonitorIdentifier> r_monid;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, 
                                          const MonitorIdentifier &monid)
{
    writeHeader(ds, r_monid, 1);
    
    SireStream::savePolyPointer(ds, monid.d);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, MonitorIdentifier &monid)
{
    VersionID v = readHeader(ds, r_monid);
    
    if (v == 1)
    {
        SireStream::loadPolyPointer(ds, monid.d);
    }
    else
        throw version_error( v, "1", r_monid, CODELOC );
        
    return ds;
}

/** Null constructor */
MonitorIdentifier::MonitorIdentifier() : MonitorID()
{}

/** Construct from the passed MonitorID */
MonitorIdentifier::MonitorIdentifier(const MonitorID &monid)
                  : MonitorID()
{
    if (monid.isA<MonitorIdentifier>())
        d = monid.asA<MonitorIdentifier>().d;
    else if (not monid.isNull())
        d.reset( monid.clone() );
}

/** Copy constructor */
MonitorIdentifier::MonitorIdentifier(const MonitorIdentifier &other)
                  : MonitorID(other), d(other.d)
{}

/** Destructor */
MonitorIdentifier::~MonitorIdentifier()
{}

/** Is this selection null? */
bool MonitorIdentifier::isNull() const
{
    return d.get() == 0;
}

/** Return a hash of this identifier */
uint MonitorIdentifier::hash() const
{
    if (d.get() == 0)
        return 0;
    else
        return d->hash();
}
            
/** Return a string representatio of this ID */
QString MonitorIdentifier::toString() const
{
    if (d.get() == 0)
        return "null";
    else
        return d->toString();
}

/** Return the base type of this ID */
const MonitorID& MonitorIdentifier::base() const
{
    if (d.get() == 0)
        return *this;
    else
        return *d;
}

/** Copy assignment operator */
MonitorIdentifier& MonitorIdentifier::operator=(const MonitorIdentifier &other)
{
    d = other.d;
    return *this;
}

/** Copy assignment operator */
MonitorIdentifier& MonitorIdentifier::operator=(const MonitorID &other)
{
    if (other.isA<MonitorIdentifier>())
        d = other.asA<MonitorIdentifier>().d;
    else if (other.isNull())
        d.reset();
    else
        d.reset(other.clone());
    
    return *this;
}

/** Comparison operator */
bool MonitorIdentifier::operator==(const SireID::ID &other) const
{
    return SireID::ID::compare<MonitorIdentifier>(*this, other);
}

/** Comparison operator */
bool MonitorIdentifier::operator==(const MonitorIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() == other.d.get();
    else
        return d == other.d or *d == *(other.d);
}

/** Comparison operator */
bool MonitorIdentifier::operator!=(const MonitorIdentifier &other) const
{
    if (d.get() == 0 or other.d.get() == 0)
        return d.get() != other.d.get();
    else
        return d != other.d and *d != *(other.d);
}

/** Comparison operator */
bool MonitorIdentifier::operator==(const MonitorID &other) const
{
    if (d.get() == 0)
        return other.isNull();
    else if (other.isA<MonitorIdentifier>())
        return this->operator==(other.asA<MonitorIdentifier>());
    else
        return d->operator==(other);
}

/** Comparison operator */
bool MonitorIdentifier::operator!=(const MonitorID &other) const
{
    if (d.get() == 0)
        return not other.isNull();
    else if (other.isA<MonitorIdentifier>())
        return this->operator!=(other.asA<MonitorIdentifier>());
    else
        return d->operator!=(other);
}

QList<MonitorName> MonitorIdentifier::map(const SystemMonitors &monitors) const
{
    if (d.get() == 0)
    {
        QList<MonitorName> names = monitors.monitorNames();
        
        if (names.isEmpty())
            throw SireSystem::missing_monitor( QObject::tr(
                "There are no monitors available in the system at all!"),
                    CODELOC );
                    
        return names;
    }
    else
        return d->map(monitors);
}

const char* MonitorIdentifier::typeName()
{
    return QMetaType::typeName( qMetaTypeId<MonitorIdentifier>() );
}

MonitorIdentifier* MonitorIdentifier::clone() const
{
    return new MonitorIdentifier(*this);
}
