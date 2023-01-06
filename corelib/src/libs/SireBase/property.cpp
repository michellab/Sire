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

#include <QMutex>

#include "property.h"
#include "propertylist.h"
#include "generalunitproperty.h"

#include <QDebug>

#include "SireError/getbacktrace.h"

#include "SireError/errors.h"
#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireBase;
using namespace SireStream;

Q_GLOBAL_STATIC_WITH_ARGS( QMutex, getGlobalMutex, ((QMutex::Recursive)) );

namespace SireBase
{

/** Return a pointer to a global mutex */
QMutex *globalLock()
{
    return getGlobalMutex();
}

} // end of namespace SireBase

///////////////
/////////////// Implementation of Property
///////////////

/** Constructor */
Property::Property() : RefCountData()
{}

/** Copy constructor */
Property::Property(const Property&) : RefCountData()
{}

/** Destructor */
Property::~Property()
{}

/** Assignment operator */
Property& Property::operator=(const Property&)
{
    return *this;
}

/** Comparison operator */
bool Property::operator==(const Property&) const
{
    return true;
}

/** Comparison operator */
bool Property::operator!=(const Property&) const
{
    return false;
}

/** Default 'toString()' function for properties - it would
    help if all properties output something more sensible */
QString Property::toString() const
{
    return QString("%1()").arg( this->what() );
}

/** Return whether or not this property holds a string (or can convert
    to a string) */
bool Property::isAString() const
{
    return false;
}

/** Return whether or not this property holds a unit (or can convert
    to a unit)
*/
bool Property::isAUnit() const
{
    return isADouble();
}

/** Return whether or not this property holds a double (or can convert
    to a double) */
bool Property::isADouble() const
{
    return false;
}

/** Return whether or not this property holds an integer (or can convert
    to an integer) */
bool Property::isAnInteger() const
{
    return false;
}

/** Return whether or not this is an array property (or can convert to an
    array property) */
bool Property::isAnArray() const
{
    return true;
}

/** Return whether or not this property holds a bool (or can convert
    to a bool) */
bool Property::isABoolean() const
{
    return false;
}

/** Return this property converted to a unit */
SireUnits::Dimension::GeneralUnit Property::asAUnit() const
{
    return SireUnits::Dimension::GeneralUnit(this->asADouble());
}

/** Return this property converted to a string. This throws an invalid
    cast if this is not possible */
QString Property::asAString() const
{
    this->throwInvalidCast("string");
    return QString();
}

/** Return this property converted to a double. This throws an invalid
    cast if this is not possible */
double Property::asADouble() const
{
    this->throwInvalidCast("double");
    return 0;
}

/** Return this property converted to an integer. This throws an invalid
    cast if this is not possible */
int Property::asAnInteger() const
{
    this->throwInvalidCast("integer");
    return 0;
}

/** Return this property converted to a bool. This throws an invalid
    cast if this is not possible */
bool Property::asABoolean() const
{
    this->throwInvalidCast("boolean");
    return false;
}

/** Return this property converted to an array property. By default, this
    automatically puts this property into a PropertyList and returns that */
PropertyList Property::asAnArray() const
{
    if (this->isA<PropertyList>())
    {
        return this->asA<PropertyList>();
    }
    else
    {
        return PropertyList(*this);
    }
}

/** Throw an invalid cast!

    \throw SireError::invalid_cast
*/
void Property::throwInvalidCast(const Property &other) const
{
    throw SireError::invalid_cast( QObject::tr(
            "Cannot cast from an object of class \"%1\" to an object "
            "of class \"%2\".")
                .arg(other.what()).arg(this->what()), CODELOC );
}

/** Throw an invalid cast!

    \throw SireError::invalid_cast
*/
void Property::throwInvalidCast(const char *typenam) const
{
    throw SireError::invalid_cast( QObject::tr(
            "Cannot cast from an object of class \"%1\" to an object "
            "of class \"%2\".")
                .arg(this->what()).arg(typenam), CODELOC );
}

static const RegisterMetaType<Property> r_propbase(MAGIC_ONLY,
                                                       "SireBase::Property");

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const Property&)
{
    writeHeader(ds, r_propbase, 0);

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, Property&)
{
    VersionID v = readHeader(ds, r_propbase);

    if (v != 0)
        throw version_error(v, "0", r_propbase, CODELOC);

    return ds;
}

///////////////
/////////////// Implementation of NullProperty
///////////////

static const RegisterMetaType<NullProperty> r_nullprop;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const NullProperty &property)
{
    writeHeader(ds, r_nullprop, 1);

    ds << static_cast<const Property&>(property);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, NullProperty &property)
{
    VersionID v = readHeader(ds, r_nullprop);

    if (v == 1)
    {
        ds >> static_cast<Property&>(property);
    }
    else
        throw version_error(v, "1", r_nullprop, CODELOC);

    return ds;
}

NullProperty::NullProperty()
             : ConcreteProperty<NullProperty,Property>()
{}

NullProperty::NullProperty(const NullProperty &other)
             : ConcreteProperty<NullProperty,Property>(other)
{}

NullProperty::~NullProperty()
{}

const char* NullProperty::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NullProperty>() );
}

/** Return the global null property */
const NullProperty& Property::null()
{
    return *create_shared_null<NullProperty>();
}

QString NullProperty::toString() const
{
    return QObject::tr("NULL");
}

///////////////
/////////////// Implementation of PropPtrBase
///////////////

static const RegisterMetaType<PropPtrBase> r_propptr( MAGIC_ONLY, NO_ROOT,
                                                      "SireBase::PropPtrBase" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const PropPtrBase &propptr)
{
    writeHeader(ds, r_propptr, 1);

    SharedDataStream sds(ds);

    sds << propptr.ptr;

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, PropPtrBase &propptr)
{
    VersionID v = readHeader(ds, r_propptr);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> propptr.ptr;
    }
    else
        throw version_error(v, "1", r_propptr, CODELOC);

    return ds;
}

/** Construct to hold a pointer to 'property' */
PropPtrBase::PropPtrBase(const Property &property)
                : ptr(property)
{}

/** Construct to hold a pointer to 'property' - this takes over
    ownership of the pointer */
PropPtrBase::PropPtrBase(Property *property)
            : ptr(property)
{}

/** Copy constructor */
PropPtrBase::PropPtrBase(const PropPtrBase &other)
                : ptr(other.ptr)
{}

/** Destructor */
PropPtrBase::~PropPtrBase()
{}

/** Copy assignment operator */
PropPtrBase& PropPtrBase::operator=(const PropPtrBase &other)
{
    ptr = other.ptr;
    return *this;
}

/** Comparison operator */
bool PropPtrBase::operator==(const PropPtrBase &other) const
{
    return ptr.constData() == other.ptr.constData() or
           ptr->equals(*(other.ptr));
}

/** Comparison operator */
bool PropPtrBase::operator!=(const PropPtrBase &other) const
{
    return ptr.constData() != other.ptr.constData() and
           not ptr->equals(*(other.ptr));
}

/** Comparison operator */
bool PropPtrBase::operator==(const Property &other) const
{
    return ptr.constData() == &other or
           ptr->equals(other);
}

/** Comparison operator */
bool PropPtrBase::operator!=(const Property &other) const
{
    return ptr.constData() != &other and
           not ptr->equals(other);
}

/** Is this a unique pointer to the object? */
bool PropPtrBase::unique() const
{
    return ptr.unique();
}

/** Detach this pointer from shared storage */
void PropPtrBase::detach()
{
    ptr.detach();
}

/** Allow automatic casting to a Property */
PropPtrBase::operator const Property&() const
{
    return *ptr;
}

/** Return a read-only reference to the object */
const Property& PropPtrBase::read() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return *ptr;
}

/** Return a writable reference to the object. This performs
    a copy-on-write test and action */
Property& PropPtrBase::edit()
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return *ptr;
}

/** Return a writable reference to the object. This performs
    a copy-on-write test and action - this is a synonym for PropPtr::edit */
Property& PropPtrBase::write()
{
    return PropPtrBase::edit();
}

bool PropPtrBase::isAString() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAString();
}

bool PropPtrBase::isAUnit() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAUnit();
}

bool PropPtrBase::isADouble() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isADouble();
}

bool PropPtrBase::isAnInteger() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAnInteger();
}

bool PropPtrBase::isABoolean() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isABoolean();
}

bool PropPtrBase::isAnArray() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAnArray();
}

QString PropPtrBase::asAString() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAString();
}

SireUnits::Dimension::GeneralUnit PropPtrBase::asAUnit() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAUnit();
}

double PropPtrBase::asADouble() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asADouble();
}

int PropPtrBase::asAnInteger() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAnInteger();
}

bool PropPtrBase::asABoolean() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asABoolean();
}

PropertyList PropPtrBase::asAnArray() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAnArray();
}

/** Throw an error as we can't cast 'got_type' into 'want_type'

    \throw SireError::invalid_cast
*/
void PropPtrBase::throwCastingError(const char *got_type, const char *want_type)
{
    throw SireError::invalid_cast( QObject::tr(
        "Cannot cast from a %1 into a %2.")
            .arg(got_type).arg(want_type), CODELOC );
}

///////////////
/////////////// Implementation of GlobalPropPtrBase
///////////////

static const RegisterMetaType<GlobalPropPtrBase> r_globalpropptr( MAGIC_ONLY, NO_ROOT,
                                                      "SireBase::GlobalPropPtrBase" );

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                        const GlobalPropPtrBase &propptr)
{
    writeHeader(ds, r_globalpropptr, 1);

    SharedDataStream sds(ds);

    sds << propptr.ptr;

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, GlobalPropPtrBase &propptr)
{
    VersionID v = readHeader(ds, r_globalpropptr);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> propptr.ptr;
    }
    else
        throw version_error(v, "1", r_globalpropptr, CODELOC);

    return ds;
}

/** Construct to hold a pointer to 'property' */
GlobalPropPtrBase::GlobalPropPtrBase(const Property &property)
                  : ptr(property)
{}

/** Construct to hold a pointer to 'property' - this takes over
    ownership of the pointer */
GlobalPropPtrBase::GlobalPropPtrBase(Property *property)
                  : ptr(property)
{}

/** Copy constructor */
GlobalPropPtrBase::GlobalPropPtrBase(const GlobalPropPtrBase &other)
                  : ptr(other.ptr)
{}

/** Destructor */
GlobalPropPtrBase::~GlobalPropPtrBase()
{}

/** Copy assignment operator */
GlobalPropPtrBase& GlobalPropPtrBase::operator=(const GlobalPropPtrBase &other)
{
    ptr = other.ptr;
    return *this;
}

/** Comparison operator */
bool GlobalPropPtrBase::operator==(const GlobalPropPtrBase &other) const
{
    return ptr.constData() == other.ptr.constData();
}

/** Comparison operator */
bool GlobalPropPtrBase::operator!=(const GlobalPropPtrBase &other) const
{
    return ptr.constData() != other.ptr.constData();
}

/** Comparison operator */
bool GlobalPropPtrBase::operator==(const Property &other) const
{
    return ptr.constData() == &other or
           ptr->equals(other);
}

/** Comparison operator */
bool GlobalPropPtrBase::operator!=(const Property &other) const
{
    return ptr.constData() != &other and
           not ptr->equals(other);
}

/** Is this a unique pointer to the object? */
bool GlobalPropPtrBase::unique() const
{
    return ptr.unique();
}

/** Allow automatic casting to a Property */
GlobalPropPtrBase::operator const Property&() const
{
    return *ptr;
}

/** Return a read-only reference to the object */
const Property& GlobalPropPtrBase::read() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return *ptr;
}
bool GlobalPropPtrBase::isAString() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAString();
}

bool GlobalPropPtrBase::isAUnit() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAUnit();
}

bool GlobalPropPtrBase::isADouble() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isADouble();
}

bool GlobalPropPtrBase::isAnInteger() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAnInteger();
}

bool GlobalPropPtrBase::isABoolean() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isABoolean();
}

bool GlobalPropPtrBase::isAnArray() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->isAnArray();
}

QString GlobalPropPtrBase::asAString() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAString();
}

SireUnits::Dimension::GeneralUnit GlobalPropPtrBase::asAUnit() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAUnit();
}

double GlobalPropPtrBase::asADouble() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asADouble();
}

int GlobalPropPtrBase::asAnInteger() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAnInteger();
}

bool GlobalPropPtrBase::asABoolean() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asABoolean();
}

PropertyList GlobalPropPtrBase::asAnArray() const
{
    BOOST_ASSERT( ptr.constData() != 0 );
    return ptr->asAnArray();
}

/** Throw an error as we can't cast 'got_type' into 'want_type'

    \throw SireError::invalid_cast
*/
void GlobalPropPtrBase::throwCastingError(const char *got_type, const char *want_type)
{
    throw SireError::invalid_cast( QObject::tr(
        "Cannot cast from a %1 into a %2.")
            .arg(got_type).arg(want_type), CODELOC );
}

////////
//////// Full instantiation of PropPtr<Property>
////////

PropPtr<Property>::PropPtr() : PropPtrBase( Property::null() )
{}

PropPtr<Property>::PropPtr(const Property &property)
                  : PropPtrBase(property)
{}

PropPtr<Property>::PropPtr(Property *property)
                  : PropPtrBase(property)
{}

PropPtr<Property>::PropPtr(const PropPtrBase &other)
                  : PropPtrBase(other)
{}

PropPtr<Property>::PropPtr(const PropPtr<Property> &other)
                  : PropPtrBase(other)
{}

PropPtr<Property>::~PropPtr()
{}

PropPtr<Property>& PropPtr<Property>::operator=(const PropPtr<Property> &other)
{
    PropPtrBase::operator=(other);
    return *this;
}

PropPtr<Property>& PropPtr<Property>::operator=(const Property &property)
{
    return this->operator=( PropPtr<Property>(property) );
}

PropPtr<Property>& PropPtr<Property>::operator=(Property *property)
{
    return this->operator=( PropPtr<Property>(property) );
}

PropPtr<Property>& PropPtr<Property>::operator=(const PropPtrBase &property)
{
    return this->operator=( PropPtr<Property>(property) );
}

const Property* PropPtr<Property>::operator->() const
{
    return &(PropPtrBase::read());
}

const Property& PropPtr<Property>::operator*() const
{
    return PropPtrBase::read();
}

const Property& PropPtr<Property>::read() const
{
    return PropPtrBase::read();
}

Property& PropPtr<Property>::edit()
{
    return PropPtrBase::edit();
}

const Property* PropPtr<Property>::data() const
{
    return &(PropPtrBase::read());
}

const Property* PropPtr<Property>::constData() const
{
    return &(PropPtrBase::read());
}

Property* PropPtr<Property>::data()
{
    return &(PropPtrBase::edit());
}

PropPtr<Property>::operator const Property&() const
{
    return PropPtrBase::read();
}

bool PropPtr<Property>::isNull() const
{
    return PropPtrBase::operator==( Property::null() );
}

PropPtr<Property> PropPtr<Property>::null()
{
    return PropPtr<Property>();
}

////////
//////// Full instantiation of GlobalPropPtr<Property>
////////

GlobalPropPtr<Property>::GlobalPropPtr() : GlobalPropPtrBase( Property::null() )
{}

GlobalPropPtr<Property>::GlobalPropPtr(const Property &property)
                        : GlobalPropPtrBase(property)
{}

GlobalPropPtr<Property>::GlobalPropPtr(Property *property)
                        : GlobalPropPtrBase(property)
{}

GlobalPropPtr<Property>::GlobalPropPtr(const GlobalPropPtrBase &other)
                        : GlobalPropPtrBase(other)
{}

GlobalPropPtr<Property>::GlobalPropPtr(const GlobalPropPtr<Property> &other)
                        : GlobalPropPtrBase(other)
{}

GlobalPropPtr<Property>::~GlobalPropPtr()
{}

GlobalPropPtr<Property>& GlobalPropPtr<Property>::operator=(
                                            const GlobalPropPtr<Property> &other)
{
    GlobalPropPtrBase::operator=(other);
    return *this;
}

GlobalPropPtr<Property>& GlobalPropPtr<Property>::operator=(const Property &property)
{
    return this->operator=( GlobalPropPtr<Property>(property) );
}

GlobalPropPtr<Property>& GlobalPropPtr<Property>::operator=(Property *property)
{
    return this->operator=( GlobalPropPtr<Property>(property) );
}

GlobalPropPtr<Property>& GlobalPropPtr<Property>::operator=(
                                                   const GlobalPropPtrBase &property)
{
    return this->operator=( GlobalPropPtr<Property>(property) );
}

const Property* GlobalPropPtr<Property>::operator->() const
{
    return &(GlobalPropPtrBase::read());
}

const Property& GlobalPropPtr<Property>::operator*() const
{
    return GlobalPropPtrBase::read();
}

const Property& GlobalPropPtr<Property>::read() const
{
    return GlobalPropPtrBase::read();
}

const Property* GlobalPropPtr<Property>::data() const
{
    return &(GlobalPropPtrBase::read());
}

const Property* GlobalPropPtr<Property>::constData() const
{
    return &(GlobalPropPtrBase::read());
}

GlobalPropPtr<Property>::operator const Property&() const
{
    return GlobalPropPtrBase::read();
}

bool GlobalPropPtr<Property>::isNull() const
{
    return GlobalPropPtrBase::operator==( Property::null() );
}

GlobalPropPtr<Property> GlobalPropPtr<Property>::null()
{
    return GlobalPropPtr<Property>();
}
