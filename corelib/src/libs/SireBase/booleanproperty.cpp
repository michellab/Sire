/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#include "SireBase/booleanproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<BooleanProperty> r_prop;

QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, const BooleanProperty &prop)
{
    writeHeader(ds, r_prop, 1);
    ds << prop.val;
    return ds;
}

QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, BooleanProperty &prop)
{
    VersionID v = readHeader(ds, r_prop);
    
    if (v == 1)
    {
        ds >> prop.val;
    }
    else
        throw version_error(v, "1", r_prop, CODELOC);
    
    return ds;
}

/** Constructor - this constructs the integer "0" */
BooleanProperty::BooleanProperty() : ConcreteProperty<BooleanProperty,Property>(), val(false)
{}

/** Construct from the passed boolean */
BooleanProperty::BooleanProperty(bool value)
               : ConcreteProperty<BooleanProperty,Property>(), val(value)
{}

/** Construct from the passed string */
BooleanProperty::BooleanProperty(const QString &value)
                : ConcreteProperty<BooleanProperty,Property>()
{
    QString v = value.toLower().simplified();
    
    if (v == "true" or v == "t" or v == "yes" or v == "y" or v == "on")
    {
        val = true;
    }
    else if (v == "false" or v == "f" or v == "no" or v == "n" or v == "off")
    {
        val = false;
    }
    else
    {
        bool ok;
        
        double dv = v.toDouble(&ok);
        
        if (ok)
        {
            val = dv;
        }
        else
        {
            throw SireError::invalid_cast( QObject::tr(
                    "Cannot convert the string '%1' to a boolean. Valid strings include "
                    "true/false, yes/no, t/f, y/n, on/off, 1/0, 1.0/0.0.")
                        .arg(value), CODELOC );
        }
    }
}

/** Construct from a VariantProperty */
BooleanProperty::BooleanProperty(const Property &other)
                : ConcreteProperty<BooleanProperty,Property>(other),
                  val( other.asABoolean() )
{}

/** Copy constructor */
BooleanProperty::BooleanProperty(const BooleanProperty &other)
               : ConcreteProperty<BooleanProperty,Property>(other), val(other.val)
{}

/** Destructor */
BooleanProperty::~BooleanProperty()
{}

const char* BooleanProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BooleanProperty>() );
}

/** Copy assignment operator */
BooleanProperty& BooleanProperty::operator=(const BooleanProperty &other)
{
    if (this != &other)
    {
        val = other.val;
    }
    
    return *this;
}

/** Comparison operator */
bool BooleanProperty::operator==(const BooleanProperty &other) const
{
    return val == other.val;
}

/** Comparison operator */
bool BooleanProperty::operator!=(const BooleanProperty &other) const
{
    return not this->operator==(other);
}

/** Return this number cast as a double */
bool BooleanProperty::value() const
{
    return val;
}

QString BooleanProperty::toString() const
{
    if (val)
        return QObject::tr("True");
    else
        return QObject::tr("False");
}

bool BooleanProperty::isAString() const
{
    return true;
}

bool BooleanProperty::isADouble() const
{
    return true;
}

bool BooleanProperty::isAnInteger() const
{
    return true;
}

bool BooleanProperty::isABoolean() const
{
    return true;
}

QString BooleanProperty::asAString() const
{
    return this->toString();
}

double BooleanProperty::asADouble() const
{
    return val;
}

int BooleanProperty::asAnInteger() const
{
    return val;
}

bool BooleanProperty::asABoolean() const
{
    return val;
}

