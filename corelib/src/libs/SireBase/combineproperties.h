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

#ifndef SIREBASE_COMBINEPROPERTIES_H
#define SIREBASE_COMBINEPROPERTIES_H

#include "property.h"
#include "propertymap.h"

#include <QList>
#include <QVector>

SIRE_BEGIN_HEADER

namespace SireBase
{
class CombineProperties;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::CombineProperties&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::CombineProperties&);

namespace SireBase
{

/** This is the base class of a property alias which provides
    a combination of existing properties. Use this to build
    properties that don't exist independently, but are
    rather built by combining together other properties.
    
    For a good example, see SireVol::CombineSpaces
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT CombineProperties : public Property
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const CombineProperties&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, CombineProperties&);

public:
    typedef QVector<PropertyName>::const_iterator const_iterator;
    typedef const_iterator iterator;

    CombineProperties();

    CombineProperties(const PropertyName &property0);

    CombineProperties(const PropertyName &property0,
                      const PropertyName &property1);
                      
    CombineProperties(const QList<PropertyName> &properties);
    CombineProperties(const QList<QString> &properties);
                      
    CombineProperties(const QVector<PropertyName> &properties);
    CombineProperties(const QVector<QString> &properties);

    CombineProperties(const CombineProperties &other);
    
    ~CombineProperties();

    static const char* typeName()
    {
        return "SireBase::CombineProperties";
    }

    const PropertyName& operator[](int i) const;

    const PropertyName& at(int i) const;

    QString toString() const;

    int nSources() const;
    int count() const;
    int size() const;

    bool isEmpty() const;

    const_iterator constBegin() const;
    const_iterator begin() const;
    
    const_iterator constEnd() const;
    const_iterator end() const;
    
    const Property& combinedProperty() const;
    
    /** Update this combined property by fetching the necessary
        properties to combine from 'properties'
        
        \throw SireBase::missing_property
        \throw SireError::invalid_cast
        \throw SireError::incompatible_error
    */
    virtual void updateFrom(const Properties &properties)=0;
    
protected:
    CombineProperties& operator=(const CombineProperties &other);
    
    bool operator==(const CombineProperties &other) const;
    bool operator!=(const CombineProperties &other) const;

    void setCombinedProperty(const Property &property);

private:
    /** The sources of the properties to be combined together */
    QVector<PropertyName> property_sources;

    /** The combined property */
    PropertyPtr combined_property;
};

}

SIRE_EXPOSE_CLASS( SireBase::CombineProperties )

SIRE_END_HEADER

#endif
