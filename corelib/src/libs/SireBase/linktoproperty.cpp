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

#include "linktoproperty.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireID;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<LinkToProperty> r_link;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const LinkToProperty &link)
{
    writeHeader(ds, r_link, 1);
    
    SharedDataStream sds(ds);
    
    sds << link.target_source << link.id_filter
        << static_cast<const Property&>(link);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, LinkToProperty &link)
{
    VersionID v = readHeader(ds, r_link);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> link.target_source >> link.id_filter
            >> static_cast<Property&>(link);
    }
    else
        throw version_error( v, "1", r_link, CODELOC );

    return ds;
}

/** Constructor */
LinkToProperty::LinkToProperty() : ConcreteProperty<LinkToProperty,Property>()
{}

/** Construct to link to the property 'source' */
LinkToProperty::LinkToProperty(const PropertyName &source)
               : ConcreteProperty<LinkToProperty,Property>(),
                 target_source(source)
{}

/** Construct to link to the property 'source' in the objects
    identified by 'filter' */
LinkToProperty::LinkToProperty(const PropertyName &source, const ID &filter)
               : ConcreteProperty<LinkToProperty,Property>(),
                 target_source(source), id_filter(filter)
{}

/** Copy constructor */
LinkToProperty::LinkToProperty(const LinkToProperty &other)
               : ConcreteProperty<LinkToProperty,Property>(other),
                 target_source(other.target_source), id_filter(other.id_filter)
{}

/** Destructor */
LinkToProperty::~LinkToProperty()
{}

/** Copy assignment operator */
LinkToProperty& LinkToProperty::operator=(const LinkToProperty &other)
{
    if (this != &other)
    {
        target_source = other.target_source;
        id_filter = other.id_filter;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool LinkToProperty::operator==(const LinkToProperty &other) const
{
    return target_source == other.target_source and
           id_filter == other.id_filter and
           Property::operator==(other);
}

/** Comparison operator */
bool LinkToProperty::operator!=(const LinkToProperty &other) const
{
    return target_source != other.target_source or
           id_filter != other.id_filter or
           Property::operator!=(other);
}

const char* LinkToProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LinkToProperty>() );
}

/** Return a string representation of this link */
QString LinkToProperty::toString() const
{
    if (target_source.hasSource())
        return QObject::tr("LinkTo( %1 )").arg(target_source.source());
     
    else if (target_source.hasValue())
        return QObject::tr("LinkTo( %1 )").arg(target_source.value().toString());

    else
        return QObject::tr("LinkTo( NULL )");
}

/** Return the target of this link */
const PropertyName& LinkToProperty::target() const
{
    return target_source;
}

/** Return any filter for this link (this is null if there is no filter) */
const SireID::ID& LinkToProperty::filter() const
{
    return id_filter.base();
}

/** Return whether or not this link is filtered */
bool LinkToProperty::isFiltered() const
{
    return not id_filter.isNull();
}
