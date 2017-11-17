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

#include "variantproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

///////////////
/////////////// Implementation of VariantProperty
///////////////

static const RegisterMetaType<VariantProperty> r_varprop;

/** Serialise to a binary data stream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const VariantProperty &varprop)
{
    writeHeader(ds, r_varprop, 1)
        << static_cast<const Property&>(varprop)
        << static_cast<const QVariant&>(varprop);
    
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, VariantProperty &varprop)
{
    VersionID v = readHeader(ds, r_varprop);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(varprop)
           >> static_cast<QVariant&>(varprop);
    }
    else
        throw version_error(v, "1", r_varprop, CODELOC);
    
    return ds;
}

/** Null constructor */
VariantProperty::VariantProperty()
                : ConcreteProperty<VariantProperty,Property>(), QVariant()
{}

/** Construct a property equal to 'value' */
VariantProperty::VariantProperty(const QVariant &value)
                : ConcreteProperty<VariantProperty,Property>(), QVariant(value)
{}

/** Construct from a 'Property' - the property must be able to
 be cast to a VariantProperty
 
 \throw SireError::invalid_cast
 */
VariantProperty::VariantProperty(const Property &property)
                : ConcreteProperty<VariantProperty,Property>(), QVariant()
{
    *this = property;
}

VariantProperty::VariantProperty(const QString &value)
                : ConcreteProperty<VariantProperty,Property>(), QVariant(value)
{}
    
VariantProperty::VariantProperty(double value)
                : ConcreteProperty<VariantProperty,Property>(), QVariant(value)
{}

/** Copy constructor */
VariantProperty::VariantProperty(const VariantProperty &other)
: ConcreteProperty<VariantProperty,Property>(other), QVariant(other)
{}

/** Destructor */
VariantProperty::~VariantProperty()
{}

/** Throw an invalid cast error */
void VariantProperty::throwInvalidCast(const QString &typname) const
{
    throw SireError::invalid_cast( QObject::tr(
                            "Cannot convert an object of type %1 to an object of type %2.")
                                  .arg(QVariant::typeName()).arg(typname), CODELOC );
}

/** Assignment operator from a QVariant */
VariantProperty& VariantProperty::operator=(const QVariant &other)
{
    QVariant::operator=(other);
    return *this;
}

/** Assignment operator from a VariantProperty */
VariantProperty& VariantProperty::operator=(const VariantProperty &other)
{
    QVariant::operator=(other);
    Property::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool VariantProperty::operator==(const VariantProperty &other) const
{
    return QVariant::operator==(other);
}

/** Comparison operator */
bool VariantProperty::operator!=(const VariantProperty &other) const
{
    return QVariant::operator!=(other);
}

const char* VariantProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<VariantProperty>() );
}

/** String operator */
QString VariantProperty::toString() const
{
    if (this->isAString())
        return this->asAString();
    else
    {
        return QString("VariantProperty( %1() )")
                    .arg(this->typeName());
    }
}

/** Return the variant property converted to a double using the QVariant conversion function */
double VariantProperty::convertToDouble() const
{
    return this->convertTo<double>();
}

/** Return the variant property converted to an integer using the QVariant conversion function */
int VariantProperty::convertToInt() const
{
    return this->convertTo<int>();
}

/** Return the variant property converted to a string using the QVariant conversion function */
QString VariantProperty::convertToString() const
{
    return this->convertTo<QString>();
}

/** Return the variant property converted to a bool using the QVariant conversion function */
bool VariantProperty::convertToBool() const
{
    return this->convertTo<bool>();
}

/** Return whether or not this can be converted to a string */
bool VariantProperty::isAString() const
{
    return this->canConvert<QString>();
}

/** Return whether or not this can be converted to a double */
bool VariantProperty::isADouble() const
{
    return this->canConvert<double>();
}

/** Return whether or not this can be converted to an integer */
bool VariantProperty::isAnInteger() const
{
    return this->canConvert<int>();
}

/** Return whether or not this can be converted to a bool */
bool VariantProperty::isABoolean() const
{
    return this->canConvert<bool>();
}

/** Return the property converted to a string */
QString VariantProperty::asAString() const
{
    return this->convertTo<QString>();
}

/** Return the property converted to a string */
double VariantProperty::asADouble() const
{
    return this->convertTo<double>();
}

/** Return the property converted to a string */
int VariantProperty::asAnInteger() const
{
    return this->convertTo<int>();
}

/** Return the property converted to a string */
bool VariantProperty::asABoolean() const
{
    return this->convertTo<bool>();
}
