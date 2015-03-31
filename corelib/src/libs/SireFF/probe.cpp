/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "probe.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

///////////
/////////// Implementation of Probe
///////////

static const RegisterMetaType<Probe> r_probe( MAGIC_ONLY, Probe::typeName() );

QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const Probe &probe)
{
    writeHeader(ds, r_probe, 1);
    
    ds << static_cast<const Property&>(probe);
    
    return ds;
}

QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, Probe &probe)
{
    VersionID v = readHeader(ds, r_probe);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(probe);
    }
    else
        throw version_error( v, "1", r_probe, CODELOC );
        
    return ds;
}

/** Constructor */
Probe::Probe() : Property()
{}

/** Copy constructor */
Probe::Probe(const Probe &other) : Property(other)
{}

/** Destructor */
Probe::~Probe()
{}

/** Copy assignment operator */
Probe& Probe::operator=(const Probe &other)
{
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool Probe::operator==(const Probe &other) const
{
    return Property::operator==(other);
}

/** Comparison operator */
bool Probe::operator!=(const Probe &other) const
{
    return Property::operator!=(other);
}

void Probe::throwCastError(const char *desired_name) const
{
    throw SireError::invalid_cast( QObject::tr(
            "Cannot convert a Probe of type %1 to a Probe of type %2.")
                .arg(this->what()).arg(desired_name), CODELOC );
}

const char* Probe::typeName()
{
    return "SireFF::Probe";
}

Q_GLOBAL_STATIC( NullProbe, nullProbe );

const Probe& Probe::null()
{
    return *(nullProbe());
}

///////////
/////////// Implementation of NullProbe
///////////

static const RegisterMetaType<NullProbe> r_nullprobe;

QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const NullProbe &nullprobe)
{
    writeHeader(ds, r_nullprobe, 1);
    
    ds << static_cast<const Probe&>(nullprobe);
    
    return ds;
}

QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, NullProbe &nullprobe)
{
    VersionID v = readHeader(ds, r_nullprobe);
    
    if (v == 1)
    {
        ds >> static_cast<Probe&>(nullprobe);
    }
    else
        throw version_error(v, "1", r_nullprobe, CODELOC);
        
    return ds;
}

/** Constructor */
NullProbe::NullProbe() : ConcreteProperty<NullProbe,Probe>()
{}

/** Copy constructor */
NullProbe::NullProbe(const NullProbe &other) : ConcreteProperty<NullProbe,Probe>(other)
{}

/** Destructor */
NullProbe::~NullProbe()
{}

/** Copy assignment operator */
NullProbe& NullProbe::operator=(const NullProbe &other)
{
    Probe::operator=(other);
    return *this;
}

/** Comparison operator */
bool NullProbe::operator==(const NullProbe &other) const
{
    return Probe::operator==(other);
}

/** Comparison operator */
bool NullProbe::operator!=(const NullProbe &other) const
{
    return Probe::operator!=(other);
}

const char* NullProbe::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullProbe>() );
}
