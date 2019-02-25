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

#include "SireCAS/expressionproperty.h"
#include "SireCAS/values.h"

#include "SireBase/numberproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireCAS;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<ExpressionProperty> r_expprop;

QDataStream &operator<<(QDataStream &ds, const ExpressionProperty &expprop)
{
    writeHeader(ds, r_expprop, 1);
    
    SharedDataStream sds(ds);
    
    sds << expprop._val << static_cast<const Property&>(expprop);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, ExpressionProperty &expprop)
{
    VersionID v = readHeader(ds, r_expprop);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> expprop._val >> static_cast<Property&>(expprop);
    }
    else
        throw version_error(v, "1", r_expprop, CODELOC);
    
    return ds;
}

namespace SireCAS
{
    PropertyPtr wrap(const ExBase &val)
    {
        return ExpressionProperty(val);
    }
    
    PropertyPtr wrap(const Expression &val)
    {
        return ExpressionProperty(val);
    }
}

/** Constructor - this constructs an empty expression */
ExpressionProperty::ExpressionProperty() : ConcreteProperty<ExpressionProperty,Property>()
{}

/** Construct from the passed expression */
ExpressionProperty::ExpressionProperty(const Expression &value)
               : ConcreteProperty<ExpressionProperty,Property>(),
                 _val(value)
{}

/** Construct from the passed expression */
ExpressionProperty::ExpressionProperty(const ExBase &value)
               : ConcreteProperty<ExpressionProperty,Property>(),
                 _val(value)
{}

/** Copy constructor */
ExpressionProperty::ExpressionProperty(const ExpressionProperty &other)
                   : ConcreteProperty<ExpressionProperty,Property>(other),
                     _val(other._val)
{}

/** Destructor */
ExpressionProperty::~ExpressionProperty()
{}

const char* ExpressionProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ExpressionProperty>() );
}

/** Copy assignment operator */
ExpressionProperty& ExpressionProperty::operator=(const ExpressionProperty &other)
{
    if (this != &other)
    {
        _val = other._val;
    }
    
    return *this;
}

/** Comparison operator */
bool ExpressionProperty::operator==(const ExpressionProperty &other) const
{
    return this->value() == other.value();
}

/** Comparison operator */
bool ExpressionProperty::operator!=(const ExpressionProperty &other) const
{
    return not this->operator==(other);
}

/** Return this number cast as an integer */
int ExpressionProperty::asAnInteger() const
{
    if (not _val.isConstant())
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert the expression '%1' to an integer")
                .arg(_val.toString()), CODELOC );

    return NumberProperty(_val.evaluate(Values())).asAnInteger();
}

/** Return this number cast as a double */
double ExpressionProperty::asADouble() const
{
    if (not _val.isConstant())
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert the expression '%1' to a double")
                .arg(_val.toString()), CODELOC );

    return _val.evaluate(Values());
}

/** Return this number cast as a boolean */
bool ExpressionProperty::asABoolean() const
{
    if (not _val.isConstant())
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert the expression '%1' to a bool")
                .arg(_val.toString()), CODELOC );

    return _val.evaluate(Values());
}

/** Return this number cast as a double */
Expression ExpressionProperty::value() const
{
    return _val;
}

QString ExpressionProperty::toString() const
{
    return _val.toString();
}

bool ExpressionProperty::isADouble() const
{
    return _val.isConstant();
}

bool ExpressionProperty::isAnInteger() const
{
    if (_val.isConstant())
    {
        return NumberProperty(_val.evaluate(Values())).isAnInteger();
    }
    else
        return false;
}

bool ExpressionProperty::isABoolean() const
{
    return _val.isConstant();
}
