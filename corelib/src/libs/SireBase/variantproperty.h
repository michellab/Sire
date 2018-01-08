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

#ifndef SIREBASE_VARIANTPROPERTY_H
#define SIREBASE_VARIANTPROPERTY_H

#include "property.h"

SIRE_BEGIN_HEADER

namespace SireBase
{
    class VariantProperty;
}

QDataStream& operator<<(QDataStream&, const SireBase::VariantProperty&);
QDataStream& operator>>(QDataStream&, SireBase::VariantProperty&);

namespace SireBase
{

/** This is a simple property that holds any value as a QVariant. This
 is designed to be used for metadata that doesn't need any tight
 checking (e.g. the author of the molecule file, the source of
 the coordinates, the 'header' lines etc.)
 
 @author Christopher Woods
 */
class SIREBASE_EXPORT VariantProperty
            : public ConcreteProperty<VariantProperty,Property>, public QVariant
{
public:
    VariantProperty();
    
    VariantProperty(const QVariant &value);
    
    VariantProperty(const Property &other);

    VariantProperty(const QString &value);
    
    VariantProperty(double value);
    
    VariantProperty(const VariantProperty &other);
    
    virtual ~VariantProperty();
    
    VariantProperty& operator=(const QVariant &value);
    VariantProperty& operator=(const VariantProperty &other);
    
    bool operator==(const VariantProperty &other) const;
    bool operator!=(const VariantProperty &other) const;
    
    static const char* typeName();
    
    QString toString() const;
    
    double convertToDouble() const;
    int convertToInt() const;
    QString convertToString() const;
    bool convertToBool() const;
    
    template<class T>
    T convertTo() const
    {
        if (not this->canConvert<T>())
            this->throwInvalidCast( QMetaType::typeName( qMetaTypeId<T>() ) );
        
        return this->value<T>();
    }
    
    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;
    
    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    
private:
    void throwInvalidCast(const QString &typname) const;
};

}

Q_DECLARE_METATYPE( SireBase::VariantProperty )

SIRE_EXPOSE_CLASS( SireBase::VariantProperty )

SIRE_END_HEADER

#endif
