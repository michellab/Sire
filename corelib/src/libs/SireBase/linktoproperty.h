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

#ifndef SIREBASE_LINKTOPROPERTY_H
#define SIREBASE_LINKTOPROPERTY_H

#include "property.h"
#include "propertymap.h"

#include "SireID/identifier.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class LinkToProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::LinkToProperty&);
QDataStream& operator>>(QDataStream&, SireBase::LinkToProperty&);

namespace SireBase
{

/** This property is actually an alias for another property.

    It allows a single property to be referenced by multiple
    different names. Also, optional identifiers can be used
    so that this link only applies to properties in specifically
    identified objects.
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT LinkToProperty : public ConcreteProperty<LinkToProperty,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const LinkToProperty&);
friend QDataStream& ::operator>>(QDataStream&, LinkToProperty&);

public:
    LinkToProperty();
    
    LinkToProperty(const PropertyName &source);
    LinkToProperty(const PropertyName &source, const SireID::ID &filter);
    
    LinkToProperty(const LinkToProperty &other);
    
    ~LinkToProperty();
    
    LinkToProperty& operator=(const LinkToProperty &other);
    
    bool operator==(const LinkToProperty &other) const;
    bool operator!=(const LinkToProperty &other) const;
    
    static const char* typeName();

    QString toString() const;

    const PropertyName& target() const;
    
    const SireID::ID& filter() const;

    bool isFiltered() const;

private:
    /** The target of this alias */
    PropertyName target_source;
    
    /** The object identity filter used to select
        only properties that match identified objects */
    SireID::Identifier id_filter;
};

}

Q_DECLARE_METATYPE( SireBase::LinkToProperty )

SIRE_EXPOSE_CLASS( SireBase::LinkToProperty )

SIRE_END_HEADER

#endif
