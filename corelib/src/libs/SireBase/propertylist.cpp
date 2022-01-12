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

#include "propertylist.h"
#include "stringproperty.h"
#include "numberproperty.h"
#include "arrayproperty.hpp"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "tostring.h"

using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<PropertyList> r_proplist;

/////////
///////// Implementation of SireBase::wrap
/////////

namespace SireBase
{
    namespace detail
    {
        int checkIndex(int i, int count)
        {
            int idx = i;

            if (i < 0)
                idx = count + i;

            if (idx < 0 or idx >= count)
                throw SireError::invalid_index( QObject::tr(
                        "Cannot access element %1. The number of elements is %2.")
                            .arg(i).arg(count), CODELOC );

            return idx;
        }
    }

    PropertyPtr wrap(const Property &value)
    {
        return PropertyPtr(value);
    }

    PropertyPtr wrap(const QList<PropertyPtr> &value)
    {
        return PropertyList(value);
    }

    PropertyPtr wrap(const QString &value)
    {
        return StringProperty(value);
    }

    PropertyPtr wrap(double value)
    {
        return NumberProperty(value);
    }

    PropertyPtr wrap(const QList<int> &values)
    {
        QVector<qint64> ivals;
        ivals.reserve(values.count());

        foreach(int value, values)
        {
            ivals.append(value);
        }

        return IntegerArrayProperty(ivals);
    }

    PropertyPtr wrap(const QList<double> &values)
    {
        return DoubleArrayProperty(values);
    }

    PropertyPtr wrap(const QVector<int> &values)
    {
        QVector<qint64> ivals;
        ivals.reserve(values.count());

        foreach(int value, values)
        {
            ivals.append(value);
        }

        return IntegerArrayProperty(ivals);
    }

    PropertyPtr wrap(const QVector<double> &values)
    {
        return DoubleArrayProperty(values);
    }

    PropertyPtr wrap(const QList<QString> &values)
    {
        return StringArrayProperty(values);
    }

    PropertyPtr wrap(const QVector<QString> &values)
    {
        return StringArrayProperty(values);
    }

    PropertyPtr wrap(const QStringList &values)
    {
        return StringArrayProperty(values);
    }
}

////////
//////// Implementation of PropertyList
////////

QDataStream &operator<<(QDataStream &ds, const PropertyList &list)
{
    writeHeader(ds, r_proplist, 1);

    SharedDataStream sds(ds);
    sds << list.l;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, PropertyList &list)
{
    VersionID v = readHeader(ds, r_proplist);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> list.l;
    }
    else
        throw version_error(v, "1", r_proplist, CODELOC);

    return ds;
}

/** Constructor */
PropertyList::PropertyList() : ConcreteProperty<PropertyList,Property>()
{}

/** Construct from the passed list */
PropertyList::PropertyList(const QList<double> &numbers)
             : ConcreteProperty<PropertyList,Property>()
{
    for (auto val : numbers)
    {
        l.append( NumberProperty(val) );
    }
}


/** Construct from the passed list */
PropertyList::PropertyList(const QStringList &strings)
             : ConcreteProperty<PropertyList,Property>()
{
    for (auto string : strings)
    {
        l.append( StringProperty(string) );
    }
}

/** Construct from the passed list */
PropertyList::PropertyList(const DoubleArrayProperty &array)
             : ConcreteProperty<PropertyList,Property>()
{
    for (auto val : array.toVector())
    {
        l.append( NumberProperty(val) );
    }
}

/** Construct from the passed list */
PropertyList::PropertyList(const IntegerArrayProperty &array)
             : ConcreteProperty<PropertyList,Property>()
{
    for (auto val : array.toVector())
    {
        l.append( NumberProperty(val) );
    }
}

/** Construct from the passed list */
PropertyList::PropertyList(const StringArrayProperty &array)
             : ConcreteProperty<PropertyList,Property>()
{
    for (auto val : array.toVector())
    {
        l.append( StringProperty(val) );
    }
}

/** Construct from the passed list */
PropertyList::PropertyList(const QList<PropertyPtr> &props)
             : ConcreteProperty<PropertyList,Property>(), l(props)
{}

/** Construct to hold the passed single value */
PropertyList::PropertyList(const Property &other)
             : ConcreteProperty<PropertyList,Property>()
{
    l.append(other);
}

/** Copy constructor */
PropertyList::PropertyList(const PropertyList &other)
             : ConcreteProperty<PropertyList,Property>(other), l(other.l)
{}

/** Destructor */
PropertyList::~PropertyList()
{}

QString PropertyList::toString() const
{
    return Sire::toString(l);
}

QList<PropertyPtr> PropertyList::array() const
{
    return l;
}

/** Copy assignment operator */
PropertyList& PropertyList::operator=(const PropertyList &other)
{
    l = other.l;
    return *this;
}

/** Comparison operator */
bool PropertyList::operator==(const PropertyList &other) const
{
    return l == other.l;
}

/** Comparison operator */
bool PropertyList::operator!=(const PropertyList &other) const
{
    return l != other.l;
}

/** Addition operator */
PropertyList PropertyList::operator+(const PropertyList &other) const
{
    return PropertyList( l + other.l );
}

/** Add the passed list onto this list */
PropertyList& PropertyList::operator+=(const Property &other)
{
    l += PropertyPtr(other);
    return *this;
}

/** Access the 'ith' property */
const Property& PropertyList::operator[](int i) const
{
    int idx = i;

    if (i < 0)
        idx = l.count() + i;

    if (idx < 0 or idx >= l.count())
        throw SireError::invalid_index( QObject::tr(
                    "Cannot access element at index %1. Number of elements in list equals %2.")
                        .arg(i).arg(l.count()), CODELOC );

    return l.at(idx).read();
}

const char* PropertyList::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PropertyList>() );
}

/** Return the number of elements in the list */
int PropertyList::count() const
{
    return l.count();
}

/** Return the number of elements in the list */
int PropertyList::size() const
{
    return l.count();
}

/** Add the passed property onto the end of this list */
void PropertyList::append(const Property &property)
{
    l.append( PropertyPtr(property) );
}

/** Append the list of properties onto the end of this list */
void PropertyList::append(const QList<PropertyPtr> &props)
{
    foreach (const PropertyPtr &prop, props)
        l.append(prop);
}

/** Return the element at index i */
const Property& PropertyList::at(int i) const
{
    return this->operator[](i);
}

/** Clear the list */
void PropertyList::clear()
{
    l.clear();
}

/** Return whether or not the list is empty */
bool PropertyList::empty() const
{
    return l.empty();
}

/** Return whether or not the list is empty */
bool PropertyList::isEmpty() const
{
    return l.isEmpty();
}

/** Insert the passed value at index 'i' */
void PropertyList::insert(int i, const Property &value)
{
    l.insert(i, PropertyPtr(value));
}

/** Return the sub-set of this list from list[pos] to list[pos+length]. If
    length is -1 then the whole rest of the list is returned */
PropertyList PropertyList::mid(int pos, int length) const
{
    return PropertyList(l.mid(pos,length));
}

/** Move an element of the list from index 'from' to index 'to' */
void PropertyList::move(int from, int to)
{
    l.move(from, to);
}

/** Pop off an element from the back of the list */
void PropertyList::pop_back()
{
    l.pop_back();
}

/** Pop off an element from the front of the list */
void PropertyList::pop_front()
{
    l.pop_front();
}

/** Prepend the passed property to the beginning of the list */
void PropertyList::prepend(const Property &value)
{
    l.prepend(value);
}

/** Push the passed value onto the back of the list */
void PropertyList::push_back(const Property &value)
{
    l.push_back(value);
}

/** Push the passed value onto the front of the list */
void PropertyList::push_front(const Property &value)
{
    l.push_front(value);
}

/** Remove the item at index 'i' */
void PropertyList::removeAt(int i)
{
    int idx = i;

    if (i < 0)
        idx = l.count() + i;

    if (idx < 0 or idx >= l.count())
        throw SireError::invalid_index( QObject::tr(
                    "Cannot remove element at index %1. Number of elements in list equals %2.")
                        .arg(i).arg(l.count()), CODELOC );

    l.removeAt(idx);
}

/** Remove the first element from the list */
void PropertyList::removeFirst()
{
    l.removeFirst();
}

/** Remove the last element from the list */
void PropertyList::removeLast()
{
    l.removeLast();
}

/** Replace the element at index 'i' with 'value' */
void PropertyList::replace(int i, const Property &value)
{
    int idx = i;

    if (i < 0)
        idx = l.count() + i;

    if (idx < 0 or idx >= l.count())
        throw SireError::invalid_index( QObject::tr(
                    "Cannot replace element at index %1. Number of elements in list equals %2.")
                        .arg(i).arg(l.count()), CODELOC );

    l.replace(idx, value);
}

/** Swap this list with 'other' */
void PropertyList::swap(PropertyList &other)
{
    QList<PropertyPtr> tmp = l;
    l = other.l;
    other.l = tmp;
}

/** Swap elements 'i' and 'j' */
void PropertyList::swap(int i, int j)
{
    int idx_i = i;

    if (i < 0)
        idx_i = l.count() + i;

    int idx_j = j;

    if (j < 0)
        idx_j = l.count() + j;

    if (idx_i < 0 or idx_j < 0 or idx_i >= l.count() or idx_j >= l.count())
        throw SireError::invalid_index( QObject::tr(
                    "Cannot swap elements %1 and %2. Number of elements in list equals %3.")
                        .arg(i).arg(j).arg(l.count()), CODELOC );

    l.swapItemsAt(idx_i, idx_j);
}

/** Take the element at index 'i' */
PropertyPtr PropertyList::takeAt(int i)
{
    int idx = i;

    if (i < 0)
        idx = l.count() + i;

    if (idx < 0 or idx >= l.count())
        throw SireError::invalid_index( QObject::tr(
                    "Cannot take element at index %1. Number of elements in list equals %2.")
                        .arg(i).arg(l.count()), CODELOC );

    return l.takeAt(idx);
}

/** Take the first element */
PropertyPtr PropertyList::takeFirst()
{
    return l.takeFirst();
}

/** Take the last element */
PropertyPtr PropertyList::takeLast()
{
    return l.takeLast();
}

/** Return this as a QList<PropertyPtr> */
QList<PropertyPtr> PropertyList::toList() const
{
    return l;
}

/** Return this as a QVector<PropertyPtr> */
QVector<PropertyPtr> PropertyList::toVector() const
{
    return l.toVector();
}

/** Return the value at index 'i' or 'default_value' if this is an invalid index */
PropertyPtr PropertyList::value(int i, const Property &default_value) const
{
    int idx = i;

    if (i < 0)
        idx = l.count() + i;

    if (idx < 0 or idx >= l.count())
        return default_value;

    return l.at(idx);
}

/** Return the value at index i, or a null property if this is
    an invalid index */
PropertyPtr PropertyList::value(int i) const
{
    return this->value(i, PropertyPtr());
}

/** Allow casting to the underlying QList<PropertyPtr> */
PropertyList::operator QList<PropertyPtr>() const
{
    return l;
}

bool PropertyList::isAString() const
{
    if (l.count() == 1)
        return l.at(0).read().isAString();
    else
        return false;
}

bool PropertyList::isADouble() const
{
    if (l.count() == 1)
        return l.at(0).read().isADouble();
    else
        return false;
}

bool PropertyList::isAnInteger() const
{
    if (l.count() == 1)
        return l.at(0).read().isAnInteger();
    else
        return false;
}

bool PropertyList::isABoolean() const
{
    if (l.count() == 1)
        return l.at(0).read().isABoolean();
    else
        return false;
}

QString PropertyList::asAString() const
{
    if (l.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert %s to a string").arg(this->toString()), CODELOC );

    return l.at(0).read().asAString();
}

double PropertyList::asADouble() const
{
    if (l.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert %s to a double").arg(this->toString()), CODELOC );

    return l.at(0).read().asADouble();
}

int PropertyList::asAnInteger() const
{
    if (l.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert %s to an integer").arg(this->toString()), CODELOC );

    return l.at(0).read().asAnInteger();
}

bool PropertyList::asABoolean() const
{
    if (l.count() != 1)
        throw SireError::invalid_cast( QObject::tr(
            "Cannot convert %s to a boolean").arg(this->toString()), CODELOC );

    return l.at(0).read().asABoolean();
}

