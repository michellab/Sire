/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "ffdetail.h"

#include "SireBase/errors.h"
#include "SireBase/propertylist.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QMutex>

using namespace SireFF;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<FFDetail> r_ff( MAGIC_ONLY, FFDetail::typeName() );

QDataStream &operator<<(QDataStream &ds, const FFDetail &ff)
{
    writeHeader(ds, r_ff, 1);

    SharedDataStream sds(ds);
    
    sds << ff.props << static_cast<const Property&>(ff);
    
    return ds;
}

QDataStream &operator>>(QDataStream &ds, FFDetail &ff)
{
    VersionID v = readHeader(ds, r_ff);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> ff.props >> static_cast<Property&>(ff);
    }
    else
        throw version_error(v, "1", r_ff, CODELOC);
    
    return ds;
}

/** Register of all forcefields */
typedef QHash<QString,PropertyPtr> FFDetailRegistry;

Q_GLOBAL_STATIC( QMutex, registryMutex );
Q_GLOBAL_STATIC( FFDetailRegistry, getRegistry );

/** Return a list of all of the forcefields that have been registered */
QStringList FFDetail::forcefields()
{
    QMutexLocker lkr( registryMutex() );
    
    return getRegistry()->keys();
}

/** Return the forcefield that has been registered with this name. This
    returns a null property if there is no property with this name */
PropertyPtr FFDetail::get(QString forcefield)
{
    QMutexLocker lkr( registryMutex() );
    
    return getRegistry()->value(forcefield);
}

/** Register the passed forcefield. This raises an error if this forcefield
    is already registered, and it differs to what is passed */
PropertyPtr FFDetail::registerForceField(const FFDetail &ff)
{
    QMutexLocker lkr( registryMutex() );
    
    auto oldff = getRegistry()->value(ff.name());
    
    if (not oldff.isNull())
    {
        if (not oldff.read().equals(ff))
            throw SireError::invalid_key( QObject::tr(
                "You cannot have two different forcefield types registered with the "
                "same name (%1). The first forcefield is %2, while the second is %3.")
                    .arg(ff.name()).arg(oldff.read().toString()).arg(ff.toString()),
                        CODELOC );
    }
    
    getRegistry()->insert(ff.name(), ff);
    
    return getRegistry()->value(ff.name());
}

/** Null constructor */
FFDetail::FFDetail() : Property()
{}

/** Construct for the forcefield called 'name' */
FFDetail::FFDetail(const QString &name) : Property()
{
    props.setProperty("name", wrap(name));
}

/** Copy constructor */
FFDetail::FFDetail(const FFDetail &other) : Property(other), props(other.props)
{}

/** Destructor */
FFDetail::~FFDetail()
{}

/** Copy assignment operator */
FFDetail& FFDetail::operator=(const FFDetail &other)
{
    props = other.props;
    Property::operator=(other);
    return *this;
}

/** Comparison operator */
bool FFDetail::operator==(const FFDetail &other) const
{
    return props == other.props;
}

/** Comparison operator */
bool FFDetail::operator!=(const FFDetail &other) const
{
    return not FFDetail::operator==(other);
}

const char* FFDetail::typeName()
{
    return "SireFF::FFDetail";
}

/** Return whether or not this is null */
bool FFDetail::isNull() const
{
    return props.isEmpty();
}

/** Return the name of the forcefield */
QString FFDetail::name() const
{
    if (isNull())
        return "NULL";
    else
        return props.property("name").asAString();
}

/** Return all of the properties of this forcefield type */
Properties FFDetail::properties() const
{
    return props;
}

/** Internal function used to set a property. Note that a property cannot
    be set twice! */
void FFDetail::setProperty(const QString &key, const Property &value)
{
    if (props.hasProperty(key))
        throw SireBase::duplicate_property( QObject::tr(
            "You cannot set a value of the property '%1' twice! Current value "
            "is %2, while the new value is %3.")
                .arg(key).arg(props.property(key).toString()).arg(value.toString()),
                    CODELOC );
    
    props.setProperty(key, value);
}

/** Internal function used to return the current value of a property.
    This raises an error if the property does not exist */
const Property& FFDetail::property(const QString &key) const
{
    if (not props.hasProperty(key))
        throw SireBase::missing_property( QObject::tr(
            "The forcefield described by '%1' does not have a property called '%2'.")
                .arg(this->toString()).arg(key), CODELOC );
    
    return props.property(key);
}

/** Assert that this forcefield is compatible with 'other' */
void FFDetail::assertCompatibleWith(const FFDetail &other) const
{
    if (not this->isCompatibleWith(other))
        throw SireError::incompatible_error( QObject::tr(
            "The two forcefields are not compatible with one another.\n%1\nversus\n%2.")
                .arg(this->toString()).arg(other.toString()), CODELOC );
}
