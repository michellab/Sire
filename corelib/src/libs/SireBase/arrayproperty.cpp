/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#include "SireBase/propertylist.h"
#include "SireBase/numberproperty.h"
#include "SireBase/stringproperty.h"
#include "SireBase/booleanproperty.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

//////////
////////// Implementation of DoubleArrayProperty
//////////

static const RegisterMetaType<DoubleArrayProperty> r_doublearray;

QDataStream &operator<<(QDataStream &ds, const DoubleArrayProperty &array)
{
    writeHeader(ds, r_doublearray, 1);

    SharedDataStream sds(ds);
    sds << array.array();

    return ds;
}

QDataStream &operator>>(QDataStream &ds, DoubleArrayProperty &array)
{
    VersionID v = readHeader(ds, r_doublearray);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<double> a;
        sds >> a;
        array = DoubleArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_doublearray, CODELOC);

    return ds;
}

DoubleArrayProperty::DoubleArrayProperty()
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >()
{}

DoubleArrayProperty::DoubleArrayProperty(const QList<double> &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(array)
{}

DoubleArrayProperty::DoubleArrayProperty(const QVector<double> &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(array)
{}

DoubleArrayProperty::DoubleArrayProperty(const StringArrayProperty &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >()
{
    for (const auto &value : array.toVector())
    {
        bool ok;
        a.append(value.toDouble(&ok));

        if (!ok)
        {
            throw SireError::invalid_arg(QObject::tr(
                "Cannot convert the string '%1' to a double!").arg(value),
                    CODELOC);
        }
    }
}

DoubleArrayProperty::DoubleArrayProperty(const IntegerArrayProperty &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(value);
    }
}

DoubleArrayProperty::DoubleArrayProperty(const PropertyList &array)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(value.asADouble());
    }
}

DoubleArrayProperty::DoubleArrayProperty(const DoubleArrayProperty &other)
                    : ConcreteProperty<DoubleArrayProperty, ArrayProperty<double> >(other)
{}

DoubleArrayProperty::~DoubleArrayProperty()
{}

const char* DoubleArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<DoubleArrayProperty>());
}

DoubleArrayProperty& DoubleArrayProperty::operator=(const DoubleArrayProperty &other)
{
    ArrayProperty<double>::operator=(other);
    return *this;
}

bool DoubleArrayProperty::operator==(const DoubleArrayProperty &other) const
{
    return ArrayProperty<double>::operator==(other);
}

bool DoubleArrayProperty::operator!=(const DoubleArrayProperty &other) const
{
    return ArrayProperty<double>::operator!=(other);
}

DoubleArrayProperty DoubleArrayProperty::operator+(const DoubleArrayProperty &other) const
{
    return DoubleArrayProperty( ArrayProperty<double>::operator+(other) );
}

DoubleArrayProperty& DoubleArrayProperty::operator+=(const DoubleArrayProperty &other)
{
    ArrayProperty<double>::operator+=(other);
    return *this;
}

bool DoubleArrayProperty::isAString() const
{
    return a.count() == 1;
}

bool DoubleArrayProperty::isADouble() const
{
    return a.count() == 1;
}

bool DoubleArrayProperty::isAnInteger() const
{
    if (a.count() == 1)
    {
        return NumberProperty(a.at(0)).isAnInteger();
    }

    return false;
}

bool DoubleArrayProperty::isABoolean() const
{
    return a.count() == 1;
}

QString DoubleArrayProperty::asAString() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a string").arg(this->toString()), CODELOC );

    return StringProperty(a.at(0)).asAString();
}

double DoubleArrayProperty::asADouble() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a double").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asADouble();
}

int DoubleArrayProperty::asAnInteger() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to an integer").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asAnInteger();
}

bool DoubleArrayProperty::asABoolean() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a boolean").arg(this->toString()), CODELOC );

    return BooleanProperty(a.at(0)).asABoolean();
}

PropertyList DoubleArrayProperty::asAnArray() const
{
    return PropertyList(*this);
}

//////////
////////// Implementation of IntegerArrayProperty
//////////

static const RegisterMetaType<IntegerArrayProperty> r_intarray;

QDataStream &operator<<(QDataStream &ds, const IntegerArrayProperty &array)
{
    writeHeader(ds, r_intarray, 1);

    SharedDataStream sds(ds);
    sds << array.array();

    return ds;
}

QDataStream &operator>>(QDataStream &ds, IntegerArrayProperty &array)
{
    VersionID v = readHeader(ds, r_intarray);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<qint64> a;
        sds >> a;
        array = IntegerArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_intarray, CODELOC);

    return ds;
}

IntegerArrayProperty::IntegerArrayProperty()
                     : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >()
{}

IntegerArrayProperty::IntegerArrayProperty(const QList<qint64> &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(array)
{}

IntegerArrayProperty::IntegerArrayProperty(const QVector<qint64> &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(array)
{}


IntegerArrayProperty::IntegerArrayProperty(const StringArrayProperty &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >()
{
    for (const auto &value : array.toVector())
    {
        bool ok;
        a.append(value.toLong(&ok));

        if (!ok)
        {
            throw SireError::invalid_arg(QObject::tr(
                "Cannot convert the string '%1' to an integer!").arg(value),
                    CODELOC);
        }
    }
}

IntegerArrayProperty::IntegerArrayProperty(const DoubleArrayProperty &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(value);
    }
}

IntegerArrayProperty::IntegerArrayProperty(const PropertyList &array)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(value.asAnInteger());
    }
}

IntegerArrayProperty::IntegerArrayProperty(const IntegerArrayProperty &other)
                    : ConcreteProperty<IntegerArrayProperty, ArrayProperty<qint64> >(other)
{}

IntegerArrayProperty::~IntegerArrayProperty()
{}

const char* IntegerArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<IntegerArrayProperty>());
}

IntegerArrayProperty& IntegerArrayProperty::operator=(const IntegerArrayProperty &other)
{
    ArrayProperty<qint64>::operator=(other);
    return *this;
}

bool IntegerArrayProperty::operator==(const IntegerArrayProperty &other) const
{
    return ArrayProperty<qint64>::operator==(other);
}

bool IntegerArrayProperty::operator!=(const IntegerArrayProperty &other) const
{
    return ArrayProperty<qint64>::operator!=(other);
}

IntegerArrayProperty IntegerArrayProperty::operator+(const IntegerArrayProperty &other) const
{
    return IntegerArrayProperty( ArrayProperty<qint64>::operator+(other) );
}

IntegerArrayProperty& IntegerArrayProperty::operator+=(const IntegerArrayProperty &other)
{
    ArrayProperty<qint64>::operator+=(other);
    return *this;
}

bool IntegerArrayProperty::isAString() const
{
    return a.count() == 1;
}

bool IntegerArrayProperty::isADouble() const
{
    return a.count() == 1;
}

bool IntegerArrayProperty::isAnInteger() const
{
    return a.count() == 1;
}

bool IntegerArrayProperty::isABoolean() const
{
    return a.count() == 1;
}

QString IntegerArrayProperty::asAString() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a string").arg(this->toString()), CODELOC );

    return StringProperty(a.at(0)).asAString();
}

double IntegerArrayProperty::asADouble() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a double").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asADouble();
}

int IntegerArrayProperty::asAnInteger() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to an integer").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asAnInteger();
}

bool IntegerArrayProperty::asABoolean() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a boolean").arg(this->toString()), CODELOC );

    return BooleanProperty(a.at(0)).asABoolean();
}

PropertyList IntegerArrayProperty::asAnArray() const
{
    return PropertyList(*this);
}

//////////
////////// Implementation of StringArrayProperty
//////////

static const RegisterMetaType<StringArrayProperty> r_stringarray;

QDataStream &operator<<(QDataStream &ds, const StringArrayProperty &array)
{
    writeHeader(ds, r_stringarray, 1);

    SharedDataStream sds(ds);
    sds << array.array();

    return ds;
}

QDataStream &operator>>(QDataStream &ds, StringArrayProperty &array)
{
    VersionID v = readHeader(ds, r_stringarray);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        QVector<QString> a;
        sds >> a;
        array = StringArrayProperty(a);
    }
    else
        throw version_error(v, "1", r_stringarray, CODELOC);

    return ds;
}

StringArrayProperty::StringArrayProperty()
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >()
{}

StringArrayProperty::StringArrayProperty(const QList<QString> &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(array)
{}

StringArrayProperty::StringArrayProperty(const QVector<QString> &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(array)
{}

StringArrayProperty::StringArrayProperty(const DoubleArrayProperty &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(QString::number(value));
    }
}

StringArrayProperty::StringArrayProperty(const IntegerArrayProperty &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(QString::number(value));
    }
}

StringArrayProperty::StringArrayProperty(const PropertyList &array)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >()
{
    for (const auto &value : array.toVector())
    {
        a.append(value.asAString());
    }
}

StringArrayProperty::StringArrayProperty(const StringArrayProperty &other)
                    : ConcreteProperty<StringArrayProperty, ArrayProperty<QString> >(other)
{}

StringArrayProperty::~StringArrayProperty()
{}

const char* StringArrayProperty::typeName()
{
    return QMetaType::typeName(qMetaTypeId<StringArrayProperty>());
}

StringArrayProperty& StringArrayProperty::operator=(const StringArrayProperty &other)
{
    ArrayProperty<QString>::operator=(other);
    return *this;
}

bool StringArrayProperty::operator==(const StringArrayProperty &other) const
{
    return ArrayProperty<QString>::operator==(other);
}

bool StringArrayProperty::operator!=(const StringArrayProperty &other) const
{
    return ArrayProperty<QString>::operator!=(other);
}

StringArrayProperty StringArrayProperty::operator+(const StringArrayProperty &other) const
{
    return StringArrayProperty( ArrayProperty<QString>::operator+(other) );
}

StringArrayProperty& StringArrayProperty::operator+=(const StringArrayProperty &other)
{
    ArrayProperty<QString>::operator+=(other);
    return *this;
}

bool StringArrayProperty::isAString() const
{
    return a.count() == 1;
}

bool StringArrayProperty::isADouble() const
{
    if (a.count() == 1)
    {
        try
        {
            return NumberProperty(a.at(0)).isADouble();
        }
        catch(...)
        {}
    }

    return false;
}

bool StringArrayProperty::isAnInteger() const
{
    if (a.count() == 1)
    {
        try
        {
            return NumberProperty(a.at(0)).isAnInteger();
        }
        catch(...)
        {}
    }

    return false;
}

bool StringArrayProperty::isABoolean() const
{
    if (a.count() == 1)
    {
        try
        {
            return BooleanProperty(a.at(0)).isABoolean();
        }
        catch(...)
        {}
    }

    return false;
}

QString StringArrayProperty::asAString() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a string").arg(this->toString()), CODELOC );

    return StringProperty(a.at(0)).asAString();
}

double StringArrayProperty::asADouble() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a double").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asADouble();
}

int StringArrayProperty::asAnInteger() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to an integer").arg(this->toString()), CODELOC );

    return NumberProperty(a.at(0)).asAnInteger();
}

bool StringArrayProperty::asABoolean() const
{
    if (a.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot cast %s to a boolean").arg(this->toString()), CODELOC );

    return BooleanProperty(a.at(0)).asABoolean();
}

PropertyList StringArrayProperty::asAnArray() const
{
    return PropertyList(*this);
}

