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

#include "combineproperties.h"
#include "properties.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireID;
using namespace SireStream;

static const RegisterMetaType<CombineProperties> r_combineprops( MAGIC_ONLY,
                                                    CombineProperties::typeName() );
                                                    
/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const CombineProperties &combineprops)
{
    writeHeader(ds, r_combineprops, 1);
    
    SharedDataStream sds(ds);
    
    sds << combineprops.property_sources << combineprops.combined_property
        << static_cast<const Property&>(combineprops);
        
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds,
                                        CombineProperties &combineprops)
{
    VersionID v = readHeader(ds, r_combineprops);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> combineprops.property_sources >> combineprops.combined_property
            >> static_cast<Property&>(combineprops);
    }
    else
        throw version_error( v, "1", r_combineprops, CODELOC );
        
    return ds;
}

/** Constructor */
CombineProperties::CombineProperties() : Property()
{}

/** Construct for just the single passed property */
CombineProperties::CombineProperties(const PropertyName &source)
                  : Property()
{
    if (not source.isNull())
    {
        property_sources = QVector<PropertyName>(1, source);
        property_sources.squeeze();
    }
}

static bool hasSource(const QString &string)
{
    return true;
}

static QString getSource(const QString &string)
{
    return string;
}

static bool hasSource(const PropertyName &name)
{
    return name.hasSource();
}

static QString getSource(const PropertyName &name)
{
    return name.source();
}

template<class T>
QVector<PropertyName> combine_properties(const T &props)
{
    QHash<QString,PropertyName> properties;
    
    for (typename T::const_iterator it = props.constBegin();
         it != props.constEnd();
         ++it)
    {
        if ( ::hasSource(*it) )
        {
            if (not properties.contains( ::getSource(*it) ))
                properties.insert( ::getSource(*it), *it );
        }
    }
    
    QVector<PropertyName> p = properties.values().toVector();
    
    if (not p.isEmpty())
        p.squeeze();
    
    return p;
}

/** Construct to combine just the passed two properties */
CombineProperties::CombineProperties(const PropertyName &property0,
                                     const PropertyName &property1)
                  : Property()
{
    QVector<PropertyName> props(2);
    props[0] = property0;
    props[1] = property1;
    
    property_sources = ::combine_properties(props);
}
                  
/** Construct to combine together the list of passed properties */
CombineProperties::CombineProperties(const QList<PropertyName> &properties)
                  : Property(), property_sources( ::combine_properties(properties) )
{}

/** Construct to combine together the list of passed properties */
CombineProperties::CombineProperties(const QList<QString> &properties)
                  : Property(), property_sources( ::combine_properties(properties) )
{}

/** Construct to combine together the list of passed properties */
CombineProperties::CombineProperties(const QVector<PropertyName> &properties)
                  : Property(), property_sources( ::combine_properties(properties) )
{}

/** Construct to combine together the list of passed properties */
CombineProperties::CombineProperties(const QVector<QString> &properties)
                  : Property(), property_sources( ::combine_properties(properties) )
{}

/** Copy constructor */
CombineProperties::CombineProperties(const CombineProperties &other)
                  : Property(other), 
                    property_sources(other.property_sources),
                    combined_property(other.combined_property)
{}

/** Destructor */
CombineProperties::~CombineProperties()
{}

/** Copy assignment operator */
CombineProperties& CombineProperties::operator=(const CombineProperties &other)
{
    if (this != &other)
    {
        property_sources = other.property_sources;
        combined_property = other.combined_property;
        Property::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool CombineProperties::operator==(const CombineProperties &other) const
{
    return this == &other or
           ( property_sources == other.property_sources and
             combined_property == other.combined_property and
             Property::operator==(other) );
}

/** Comparison operator */
bool CombineProperties::operator!=(const CombineProperties &other) const
{
    return not this->operator==(other);
}

/** Return the ith property source 

    \throw SireID::invalid_index
*/
const PropertyName& CombineProperties::operator[](int i) const
{
    return property_sources.constData()[ Index(i).map(property_sources.count()) ];
}

/** Return the ith property source 

    \throw SireID::invalid_index
*/
const PropertyName& CombineProperties::at(int i) const
{
    return this->operator[](i);
}

/** Return the number of properties that are combined together
    to form this property */
int CombineProperties::nSources() const
{
    return property_sources.count();
}

/** Return the number of properties that are combined together
    to form this property */
int CombineProperties::count() const
{
    return this->nSources();
}

/** Return the number of properties that are combined together
    to form this property */
int CombineProperties::size() const
{
    return this->nSources();
}

/** Return whether or not this is empty (has no properties!) */
bool CombineProperties::isEmpty() const
{
    return property_sources.isEmpty();
}

/** Return a string representation of this combination */
QString CombineProperties::toString() const
{
    return QString( "%1( nSources() == %2 )" )
                .arg(this->what()).arg(this->nSources());
}

CombineProperties::const_iterator CombineProperties::constBegin() const
{
    return property_sources.constBegin();
}

CombineProperties::const_iterator CombineProperties::begin() const
{
    return property_sources.begin();
}

CombineProperties::const_iterator CombineProperties::constEnd() const
{
    return property_sources.constEnd();
}

CombineProperties::const_iterator CombineProperties::end() const
{
    return property_sources.end();
}

/** Return the combined property. This will be null if this property
    has not been updated, or if there are no properties to combine */
const Property& CombineProperties::combinedProperty() const
{
    return combined_property;
}

/** Internal function used to set the combined property */
void CombineProperties::setCombinedProperty(const Property &property)
{
    combined_property = property;
}
