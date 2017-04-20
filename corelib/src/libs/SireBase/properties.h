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

#ifndef SIREBASE_PROPERTIES_H
#define SIREBASE_PROPERTIES_H

#include "property.h"
#include "propertymap.h"

#include "shareddatapointer.hpp"

SIRE_BEGIN_HEADER

namespace SireBase
{
class Properties;
}

QDataStream& operator<<(QDataStream&, const SireBase::Properties&);
QDataStream& operator>>(QDataStream&, SireBase::Properties&);

XMLStream& operator<<(XMLStream&, const SireBase::Properties&);
XMLStream& operator>>(XMLStream&, SireBase::Properties&);

QTextStream& operator<<(QTextStream&, const SireBase::Properties&);

namespace SireBase
{

namespace detail
{
class PropertiesData;
}

/** This class holds a collection of properties, indexed by name.
    Each property comes complete with a set of metadata.
    The metadata is actually another Properties object,
    and indeed Properties can itself be a Property,
    so allowing Properties to be nested indefinitely.

    @author Christopher Woods
*/
class SIREBASE_EXPORT Properties : public ConcreteProperty<Properties,Property>
{

friend QDataStream& ::operator<<(QDataStream&, const Properties&);
friend QDataStream& ::operator>>(QDataStream&, Properties&);

friend XMLStream& ::operator<<(XMLStream&, const Properties&);
friend XMLStream& ::operator>>(XMLStream&, Properties&);

friend class detail::PropertiesData; // so can call private constructor

public:
    typedef QHash<QString,PropertyPtr>::const_iterator const_iterator;
    typedef QHash<QString,PropertyPtr>::const_iterator iterator;

    Properties();

    Properties(const Properties &other);

    ~Properties();

    Properties& operator=(const Properties &other);

    bool operator==(const Properties &other) const;
    bool operator!=(const Properties &other) const;

    static const char* typeName()
    {
        return "SireBase::Properties";
    }

    const Property& operator[](const PropertyName &key) const;

    bool isEmpty() const;

    QString toString() const;

    QStringList propertyKeys() const;
    
    QStringList metadataKeys() const;
    QStringList metadataKeys(const PropertyName &key) const;

    template<class T>
    QStringList propertyKeysOfType() const;

    template<class T>
    QStringList metadataKeysOfType() const;
    
    template<class T>
    QStringList metadataKeysOfType(const PropertyName &key) const;

    const_iterator begin() const;
    const_iterator constBegin() const;
    
    const_iterator find(const QString &key) const;
    const_iterator constFind(const QString &key) const;
    
    const_iterator end() const;
    const_iterator constEnd() const;

    int count() const;
    int size() const;
    int nProperties() const;

    const Property& property(const PropertyName &key) const;
    const Property& property(const PropertyName &key,
                             const Property &default_value) const;
    
    const Property& metadata(const PropertyName &metakey) const;
    const Property& metadata(const PropertyName &metakey,
                             const Property &default_value) const;
    
    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey) const;
    const Property& metadata(const PropertyName &key,
                             const PropertyName &metakey,
                             const Property &default_value) const;

    const Properties& allMetadata() const;
    const Properties& allMetadata(const PropertyName &key) const;

    void setProperty(const QString &key, const Property &value);

    void setProperty(const QString &key, const Property &value,
                     bool clear_metadata);
    
    void setMetadata(const QString &metakey, const Property &value);
    void setMetadata(const QString &key,
                     const QString &metakey, const Property &value);

    void removeProperty(const QString &key);

    void removeMetadata(const QString &metakey);
    void removeAllMetadata();
    
    void removeMetadata(const QString &key, const QString &metakey);
    void removeAllMetadata(const QString &key);

    void clear();

    bool hasProperty(const PropertyName &key) const;
    
    bool hasMetadata(const PropertyName &metakey) const;
    
    bool hasMetadata(const PropertyName &key, 
                     const PropertyName &metakey) const;
                     
    template<class T>
    bool hasPropertyOfType(const PropertyName &key) const;
                     
    template<class T>
    bool hasMetadataOfType(const PropertyName &metakey) const;
    
    template<class T>
    bool hasMetadataOfType(const PropertyName &key,
                           const PropertyName &metakey) const;

    const char* propertyType(const PropertyName &key) const;
    
    const char* metadataType(const PropertyName &metakey) const;
    const char* metadataType(const PropertyName &key,
                             const PropertyName &metakey) const;
    
    void assertContainsProperty(const PropertyName &key) const;
    
    void assertContainsMetadata(const PropertyName &metakey) const;
    void assertContainsMetadata(const PropertyName &key, 
                                const PropertyName &metakey) const;

private:
    Properties(bool);

    /** Implicitly shared pointer to the data of this object */
    SharedDataPointer<detail::PropertiesData> d;
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Return the keys of all properties that are of type T */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Properties::propertyKeysOfType() const
{
    QStringList keys;

    for (Properties::const_iterator it = this->constBegin();
         it != this->constEnd();
         ++it)
    {
        if (it.value()->isA<T>())
            keys.append(it.key());
    }
    
    return keys;
}

/** Return the metakeys of all metadata of type T */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Properties::metadataKeysOfType() const
{
    QStringList metakeys;
    
    const Properties &metadata = this->allMetadata();
    
    for (Properties::const_iterator it = metadata.constBegin();
         it != metadata.constEnd();
         ++it)
    {
        if (it.value()->isA<T>())
            metakeys.append(it.key());
    }
    
    return metakeys;
}

/** Return the metakeys of all metadata of type T for the property
    with key 'key'
    
    \throw SireBase::missing_property
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
QStringList Properties::metadataKeysOfType(const PropertyName &key) const
{
    QStringList metakeys;
    const Properties &metadata = this->allMetadata(key);
    
    for (Properties::const_iterator it = metadata.constBegin();
         it != metadata.constEnd();
         ++it)
    {
        if (it.value()->isA<T>())
            metakeys.append(it.key());
    }
    
    return metakeys;
}

/** Return whether or not this molecule has a property called 'key'
    that is of type 'T' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Properties::hasPropertyOfType(const PropertyName &key) const
{
    return this->hasProperty(key) and 
           this->property(key).isA<T>();
}
                 
/** Return whether or not this molecule has some metadata at metakey
    'metakey' that is of type 'T' */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Properties::hasMetadataOfType(const PropertyName &metakey) const
{
    return this->hasMetadata(metakey) and
           this->metadata(metakey).isA<T>();
}

/** Return whether or not the property at key 'key' has some metadata
    at metakey 'metakey' that is of type 'T'
    
    \throw SireBase::missing_property
*/
template<class T>
SIRE_OUTOFLINE_TEMPLATE
bool Properties::hasMetadataOfType(const PropertyName &key,
                                   const PropertyName &metakey) const
{
    return this->hasMetadata(key, metakey) and
           this->metadata(key, metakey).isA<T>();
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireBase::Properties);

SIRE_EXPOSE_CLASS( SireBase::Properties )

SIRE_END_HEADER

#endif
