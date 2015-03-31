/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include "propertymap.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireBase;
using namespace SireStream;

////////////
//////////// Implementation of PropertyName
////////////

static const RegisterMetaType<PropertyName> r_propname(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, 
                                      const PropertyName &propname)
{
    writeHeader(ds, r_propname, 1);
    
    SharedDataStream sds(ds);
    
    sds << propname.src << propname.val << propname.value_is_default;
    
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, PropertyName &propname)
{
    VersionID v = readHeader(ds, r_propname);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> propname.src >> propname.val >> propname.value_is_default;
        
        if (propname.src.isEmpty())
        	propname.src = QString::null;
    }
    else
        throw version_error(v, "1", r_propname, CODELOC);

    return ds;
}

/** Null constructor */
PropertyName::PropertyName() : value_is_default(false)
{}

/** Construct a PropertyName that searches for the
    property using the source 'source' */
PropertyName::PropertyName(const char *source)
             : src(source), value_is_default(false)
{}

/** Construct a PropertyName that searches for the 
    property using the source 'source' */
PropertyName::PropertyName(const QString &source)
             : src(source), value_is_default(false)
{}

/** Construct a PropertyName that uses the supplied
    value, rather than searching for the property */
PropertyName::PropertyName(const Property &value)
             : val(value)
{}

/** Construct a PropertyName that searches for the property
    using the source 'source', but only if that source is
    specifically provided - otherwise the supplied default
    value of the property is used instead */
PropertyName::PropertyName(const QString &source, 
                           const Property &default_value)
             : src(source), val(default_value), value_is_default(true)
{
    BOOST_ASSERT(not source.isEmpty());
}

/** Copy constructor */
PropertyName::PropertyName(const PropertyName &other)
             : src(other.src), val(other.val), value_is_default(other.value_is_default)
{}

/** Destructor */
PropertyName::~PropertyName()
{}

/** Copy assignment operator */
PropertyName& PropertyName::operator=(const PropertyName &other)
{
    src = other.src;
    val = other.val;
    value_is_default = other.value_is_default;
    
    return *this;
}

/** Comparison operator */
bool PropertyName::operator==(const PropertyName &other) const
{
    return src == other.src and val == other.val and
           value_is_default == other.value_is_default;
}

/** Comparison operator */
bool PropertyName::operator!=(const PropertyName &other) const
{
    return src != other.src or val != other.val or 
           value_is_default != other.value_is_default;
}

const char* PropertyName::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PropertyName>() );
}

/** Return a PropertyName that says that this property is not set */
PropertyName PropertyName::none()
{
    return PropertyName();
}

/** Return whether or not the source has been set */
bool PropertyName::hasSource() const
{
    return not src.isEmpty();
}

/** Return whether or not the value has been set */
bool PropertyName::hasValue() const
{
    return not val.isNull();
}

/** Return whether or not this has a default value */
bool PropertyName::hasDefaultValue() const
{
    return value_is_default;
}

/** Return whether this property is null */
bool PropertyName::isNull() const
{
    return src.isEmpty() and val.isNull();
}

/** Return the source of the property - this is only valid
    if .hasSource() is true */
const QString& PropertyName::source() const
{
    return src;
}

/** Return the value of the property - this is only valid
    if .hasValue() is true */
const Property& PropertyName::value() const
{
    return val;
}

/** Return a string representation of this propertyname */
QString PropertyName::toString() const
{
    if (this->hasSource())
    {
        if (value_is_default)
            return QString("%1 {default: %2}").arg(src).arg(val->what());
        else
            return src;
    }
    else if (this->hasValue())
        return val->what();
    else
        return "NULL";
}

////////////
//////////// Implementation of PropertyMap
////////////

static const RegisterMetaType<PropertyMap> r_propmap(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, 
                                      const PropertyMap &propmap)
{
    writeHeader(ds, r_propmap, 1);
    
    SharedDataStream sds(ds);
    sds << propmap.propmap;
    
    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, 
                                      PropertyMap &propmap)
{
    VersionID v = readHeader(ds, r_propmap);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> propmap.propmap;
    }
    else
        throw version_error(v, "1", r_propmap, CODELOC);

    return ds;
}

/** Null constructor */
PropertyMap::PropertyMap()
{}

/** Construct a map that holds just a single PropertyName */
PropertyMap::PropertyMap(const QString &property, const PropertyName &propname)
{
    propmap.insert(property, propname);
}

/** Construct a map that holds lots of PropertyNames */
PropertyMap::PropertyMap(const QHash<QString,PropertyName> &propnames)
            : propmap(propnames)
{}

/** Copy constructor */
PropertyMap::PropertyMap(const PropertyMap &other)
            : propmap(other.propmap)
{}

/** Destructor */
PropertyMap::~PropertyMap()
{}

/** Copy assignment operator */
PropertyMap& PropertyMap::operator=(const PropertyMap &other)
{
    propmap = other.propmap;
    return *this;
}

/** Add two PropertyMaps together. This copies the properties of
    other and adds them to this, replacing any in this that 
    have the same name */
PropertyMap PropertyMap::operator+(const PropertyMap &other) const
{
    PropertyMap ret(*this);
    
    for (QHash<QString,PropertyName>::const_iterator it = other.propmap.begin();
         it != other.propmap.end();
         ++it)
    {
        ret.propmap.insert(it.key(), it.value());
    }
    
    return ret;
}

/** Comparison operator */
bool PropertyMap::operator==(const PropertyMap &other) const
{
    return propmap == other.propmap;
}

/** Comparison operator */
bool PropertyMap::operator!=(const PropertyMap &other) const
{
    return propmap != other.propmap;
}

/** Map the property called 'name' to the source or value of 
    that property */
PropertyName PropertyMap::operator[](const QString &name) const
{
    QHash<QString,PropertyName>::const_iterator
                                    it = propmap.constFind(name);
                                    
    if (it == propmap.constEnd())
    {
        return PropertyName(name);
    }
    else
    {
        return it.value();
    }
}

/** Map the property called 'name' to the source or value of 
    that property */
PropertyName PropertyMap::operator[](const char *name) const
{
    return this->operator[](QLatin1String(name));
}

/** Map the property in 'propname' to the source or value of
    that property. */
PropertyName PropertyMap::operator[](const PropertyName &propname) const
{
    if (propname.hasSource())
    {
        QHash<QString,PropertyName>::const_iterator 
                                        it = propmap.constFind(propname.source());
                                        
        if (it != propmap.constEnd())
            return it.value();
        else
            return propname;
    }
    else
        return propname;
}

const char* PropertyMap::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PropertyMap>() );
}

/** Return whether or not this map is default - if it is,
    then it doesn't specify any properties */
bool PropertyMap::isDefault() const
{
    return propmap.isEmpty();
}

/** Return whether or not this map specifies the source or value
    of the property called 'name' */
bool PropertyMap::specified(const QString &name) const
{
    return propmap.contains(name);
}

/** Return whether or not this map specifies the source or value
    of the property called 'name' */
bool PropertyMap::specified(const char *name) const
{
    return propmap.contains(name);
}

/** Return whether or not this map specifies the source or value
    of the property called 'name' */
bool PropertyMap::specified(const PropertyName &propname) const
{
    if (propname.hasSource())
        return propmap.contains(propname.source());
    else
        return propname.hasValue();
}

/** Set the property called 'name' to have the source or value
    in 'source'. This replaces any existing source or value
    for any existing property of this name in this map */
void PropertyMap::set(const QString &name, const PropertyName &source)
{
    if (not source.hasValue())
    {
        if (source.source() == name)
        {
            if (propmap.contains(name))
                propmap.remove(name);
                
            return;
        }
    }

    propmap.insert(name, source);
}

/** Return a string representation of this PropertyMap */
QString PropertyMap::toString() const
{
    QStringList items;
    
    for (QHash<QString,PropertyName>::const_iterator it = propmap.constBegin();
         it != propmap.constEnd();
         ++it)
    {
        items.append( QString("%1 == %2").arg( it.key(), it.value().toString() ) );
    }
    
    return QString("[ %1 ]").arg( items.join(", ") );
}
