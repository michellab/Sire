/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREBASE_PROPERTYMAP_H
#define SIREBASE_PROPERTYMAP_H

#include <QHash>
#include <QString>
#include <QList>

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
class PropertyName;
class PropertyMap;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::PropertyName&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::PropertyName&);

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::PropertyMap&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::PropertyMap&);

namespace SireBase
{

class Properties;

/** This class is used to store the registered name of the
    property (used as the offial name of the property
    by the code in the program, e.g. "coordinates" is the
    official name of the coordinates property, "charges"
    is the official name of the charges), the name for
    the property assigned by the user (so the user can
    say that the "charges" property is actually in the
    property called "my charges") and an overridden
    value for the property that is used instead of
    the value found in the object being queried
    (so the user can say to use a specific value of
     a property)

    This class is not used directly by the code, but
    is instead used as part of the Property::set( ) function,
    so that the user can write;

    cljff.add( mol, Property::set("charges","chgs") +
                    Property::set("ljs","ljparams") );

    The PropertyMap/PropertyName classes provide a kwargs
    like interface for the C++ classes - indeed the python
    wrappers should allow code to be written like;

    cljff.add( mol, {"charges" : "chgs", "ljs" : "ljparams"} )

    or

    cljff.add( mol, charges=="charges", ljs=="ljparams" )

    @author Christopher Woods
*/
class SIREBASE_EXPORT PropertyName
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const PropertyName&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, PropertyName&);

public:
    PropertyName();
    PropertyName(const char* source);

    PropertyName(const QString &source);
    PropertyName(const QString &source, const Property &default_value);

    PropertyName(const Property &value);

    PropertyName(const PropertyName &other);

    ~PropertyName();

    static const char* typeName();

    const char* what() const
    {
        return PropertyName::typeName();
    }

    PropertyName& operator=(const PropertyName &other);

    bool operator==(const PropertyName &other) const;
    bool operator!=(const PropertyName &other) const;

    bool hasSource() const;
    bool hasValue() const;

    bool hasDefaultValue() const;

    bool isNull() const;

    const QString& source() const;
    const Property& value() const;

    QString toString() const;

    static PropertyName none();

private:
    /** The name to use to find the property in the
        Properties container */
    QString src;

    /** The supplied or default value of the property */
    PropertyPtr val;

    /** Is the supplied value a default value? */
    bool value_is_default;
};

/** This is the class that holds the collection of user-supplied
    optional properties and their locations to functions.

    This class allows the following code to be written;

    cljff.add( mol, Property::set("charges","chgs") +
                    Property::set("ljs","ljparams") );

    The PropertyMap/PropertyName classes provide a kwargs
    like interface for the C++ classes - indeed the python
    wrappers should allow code to be written like;

    cljff.add( mol, {"charges" : "chgs", "ljs" : "ljparams"} )

    or

    cljff.add( mol, charges="charges", ljs="ljparams" )

    @author Christopher Woods
*/
class SIREBASE_EXPORT PropertyMap
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const PropertyMap&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, PropertyMap&);

public:
    PropertyMap();
    PropertyMap(const QString &property, const PropertyName &propname);

    PropertyMap(const QHash<QString,PropertyName> &propnames);

    PropertyMap(const PropertyMap &other);

    ~PropertyMap();

    static const char* typeName();

    const char* what() const
    {
        return PropertyMap::typeName();
    }

    PropertyMap& operator=(const PropertyMap &other);

    PropertyMap operator+(const PropertyMap &other) const;

    bool operator==(const PropertyMap &other) const;
    bool operator!=(const PropertyMap &other) const;

    PropertyName operator[](const char *name) const;
    PropertyName operator[](const QString &name) const;
    PropertyName operator[](const PropertyName &name) const;

    bool isDefault() const;

    bool specified(const char *name) const;
    bool specified(const QString &name) const;
    bool specified(const PropertyName &name) const;

    void set(const QString &name, const PropertyName &source);

    PropertyMap merge(const PropertyMap &other) const;

    QString toString() const;

    const QHash<QString,PropertyName> toDict() const;

private:
    /** Hash indexing all of the PropertyNames */
    QHash<QString,PropertyName> propmap;
};

}

Q_DECLARE_METATYPE(SireBase::PropertyName);
Q_DECLARE_METATYPE(SireBase::PropertyMap);

SIRE_EXPOSE_CLASS( SireBase::PropertyMap )
SIRE_EXPOSE_CLASS( SireBase::PropertyName )

SIRE_END_HEADER

#endif
