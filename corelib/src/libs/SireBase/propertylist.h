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

#ifndef SIREBASE_PROPERTYLIST_H
#define SIREBASE_PROPERTYLIST_H

#include "SireBase/property.h"
#include "SireBase/arrayproperty.hpp"

SIRE_BEGIN_HEADER

namespace SireBase
{
class PropertyList;
class DoubleArrayProperty;
class IntegerArrayProperty;
class StringArrayProperty;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::PropertyList&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::PropertyList&);

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::DoubleArrayProperty&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::DoubleArrayProperty&);

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::IntegerArrayProperty&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::IntegerArrayProperty&);

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::StringArrayProperty&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::StringArrayProperty&);

namespace SireBase
{

/** This class provides a simple list of properties. This is useful,
    e.g. if the user wants to store a list of values in a molecule.

    This class is designed to be used with StringProperty, NumberProperty
    and VectorProperty (from SireMaths) to make it easy to attach
    arbitrary data to any Sire object, e.g. storing a fixed center
    point in a molecule, saving a list of velocities with a molecule,
    etc.

    Note also that Properties already provides a Property dictionary.

    @author Christopher Woods
*/
class SIREBASE_EXPORT PropertyList : public ConcreteProperty<PropertyList,Property>
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const PropertyList&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, PropertyList&);

public:
    PropertyList();
    PropertyList(const QList<double> &numbers);
    PropertyList(const QStringList &strings);
    PropertyList(const DoubleArrayProperty &array);
    PropertyList(const IntegerArrayProperty &array);
    PropertyList(const StringArrayProperty &array);
    PropertyList(const QList<PropertyPtr> &props);
    PropertyList(const Property &other);
    PropertyList(const PropertyList &other);

    ~PropertyList();

    PropertyList& operator=(const PropertyList &other);
    bool operator==(const PropertyList &other) const;
    bool operator!=(const PropertyList &other) const;

    PropertyList operator+(const PropertyList &other) const;
    PropertyList& operator+=(const Property &other);
    const Property& operator[](int i) const;

    static const char* typeName();

    QString toString() const;

    QList<PropertyPtr> array() const;

    int count() const;
    int size() const;

    void append(const Property &property);
    void append(const QList<PropertyPtr> &props);

    const Property& at(int i) const;

    void clear();

    bool empty() const;
    bool isEmpty() const;

    void insert(int i, const Property &value);

    PropertyList mid(int pos, int length=-1) const;

    void move(int from, int to);

    void pop_back();
    void pop_front();
    void prepend(const Property &value);
    void push_back(const Property &value);
    void push_front(const Property &value);

    void removeAt(int i);
    void removeFirst();
    void removeLast();

    void replace(int i, const Property &value);

    void swap(PropertyList &other);

    void swap(int i, int j);

    PropertyPtr takeAt(int i);
    PropertyPtr takeFirst();
    PropertyPtr takeLast();

    QList<PropertyPtr> toList() const;
    QVector<PropertyPtr> toVector() const;

    PropertyPtr value(int i) const;
    PropertyPtr value(int i, const Property &default_value) const;

    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;

    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;

    operator QList<PropertyPtr>() const;

private:
    /** The actual list */
    QList<PropertyPtr> l;
};

class SIREBASE_EXPORT DoubleArrayProperty
        : public ConcreteProperty<DoubleArrayProperty,ArrayProperty<double> >
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const DoubleArrayProperty&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, DoubleArrayProperty&);

public:
    DoubleArrayProperty();
    DoubleArrayProperty(const QList<double> &array);
    DoubleArrayProperty(const QVector<double> &array);
    DoubleArrayProperty(const DoubleArrayProperty &other);

    ~DoubleArrayProperty();

    static const char* typeName();

    DoubleArrayProperty& operator=(const DoubleArrayProperty &other);

    bool operator==(const DoubleArrayProperty &other) const;
    bool operator!=(const DoubleArrayProperty &other) const;

    DoubleArrayProperty operator+(const DoubleArrayProperty &other) const;
    DoubleArrayProperty& operator+=(const DoubleArrayProperty &other);

    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;

    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    PropertyList asAnArray() const;
};

class SIREBASE_EXPORT IntegerArrayProperty
        : public ConcreteProperty<IntegerArrayProperty,ArrayProperty<qint64> >
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const IntegerArrayProperty&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, IntegerArrayProperty&);

public:
    IntegerArrayProperty();
    IntegerArrayProperty(const QList<qint64> &array);
    IntegerArrayProperty(const QVector<qint64> &array);
    IntegerArrayProperty(const IntegerArrayProperty &other);

    ~IntegerArrayProperty();

    static const char* typeName();

    IntegerArrayProperty& operator=(const IntegerArrayProperty &other);

    bool operator==(const IntegerArrayProperty &other) const;
    bool operator!=(const IntegerArrayProperty &other) const;

    IntegerArrayProperty operator+(const IntegerArrayProperty &other) const;
    IntegerArrayProperty& operator+=(const IntegerArrayProperty &other);

    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;

    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    PropertyList asAnArray() const;
};

class SIREBASE_EXPORT StringArrayProperty
        : public ConcreteProperty<StringArrayProperty,ArrayProperty<QString> >
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const StringArrayProperty&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, StringArrayProperty&);

public:
    StringArrayProperty();
    StringArrayProperty(const QList<QString> &array);
    StringArrayProperty(const QVector<QString> &array);
    StringArrayProperty(const StringArrayProperty &other);

    ~StringArrayProperty();

    static const char* typeName();

    StringArrayProperty& operator=(const StringArrayProperty &other);

    bool operator==(const StringArrayProperty &other) const;
    bool operator!=(const StringArrayProperty &other) const;

    StringArrayProperty operator+(const StringArrayProperty &other) const;
    StringArrayProperty& operator+=(const StringArrayProperty &other);

    bool isAString() const;
    bool isADouble() const;
    bool isAnInteger() const;
    bool isABoolean() const;

    QString asAString() const;
    double asADouble() const;
    int asAnInteger() const;
    bool asABoolean() const;
    PropertyList asAnArray() const;
};

SIREBASE_EXPORT PropertyPtr wrap(const Property &value);
SIREBASE_EXPORT PropertyPtr wrap(const QString &value);
SIREBASE_EXPORT PropertyPtr wrap(double value);
SIREBASE_EXPORT PropertyPtr wrap(const QVariant &value);
SIREBASE_EXPORT PropertyPtr wrap(const char *value);

SIREBASE_EXPORT PropertyPtr wrap(const QList<PropertyPtr> &value);
SIREBASE_EXPORT PropertyPtr wrap(const QList<double> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QList<int> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QVector<double> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QVector<int> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QList<QString> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QVector<QString> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QStringList &values);
SIREBASE_EXPORT PropertyPtr wrap(const QList<QVariant> &values);
SIREBASE_EXPORT PropertyPtr wrap(const QVector<QVariant> &values);

}

SIRE_EXPOSE_FUNCTION( SireBase::wrap )

Q_DECLARE_METATYPE( SireBase::DoubleArrayProperty )
Q_DECLARE_METATYPE( SireBase::IntegerArrayProperty )
Q_DECLARE_METATYPE( SireBase::StringArrayProperty )
Q_DECLARE_METATYPE( SireBase::PropertyList )

SIRE_EXPOSE_CLASS( SireBase::DoubleArrayProperty )
SIRE_EXPOSE_CLASS( SireBase::IntegerArrayProperty )
SIRE_EXPOSE_CLASS( SireBase::StringArrayProperty )
SIRE_EXPOSE_CLASS( SireBase::PropertyList )

SIRE_END_HEADER

#endif

