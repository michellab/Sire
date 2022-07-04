/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#include "bondtype.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<BondType> r_bondtype;

QDataStream &operator<<(QDataStream &ds, const BondType &b)
{
    writeHeader(ds, r_bondtype, 1);

    ds << b.bond_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BondType &b)
{
    VersionID v = readHeader(ds, r_bondtype);

    if (v == 1)
    {
        ds >> b.bond_type;
    }
    else
        throw version_error(v, "1", r_bondtype, CODELOC);

    return ds;
}

/** Constructor (default is an undefined bond) */
BondType::BondType() : ConcreteProperty<BondType,Property>(), bond_type(0)
{}

/** Construct from the passed string */
BondType::BondType(const QString &str)
         : ConcreteProperty<BondType,Property>()
{
    auto s = str.trimmed().toLower();

    if (s == "single")
        this->bond_type = 1;
    else if (s == "double")
        this->bond_type = 2;
    else if (s == "triple")
        this->bond_type = 3;
    else if (s == "aromatic")
        this->bond_type = 4;
    else if (s == "undefined")
        this->bond_type = 0;
    else
        throw SireError::invalid_arg(QObject::tr(
            "Cannot interpret bond type '%1'. Should be one of "
            "'single', 'double', 'triple', 'aromatic' or 'undefined'.")
              .arg(str), CODELOC);
}

/** Construct from the the passed number */
BondType::BondType(int value)
         : ConcreteProperty<BondType,Property>()
{
    if (value < 0 or value > 4)
        throw SireError::invalid_arg(QObject::tr(
            "Invalid bond type '%1'. Should be an integer between "
            "0 and 4.").arg(value), CODELOC);

    this->bond_type = value;
}

/** Copy constructor */
BondType::BondType(const BondType &other)
         : ConcreteProperty<BondType,Property>(other),
           bond_type(other.bond_type)
{}

/** Destructor */
BondType::~BondType()
{}

/** Copy assignment operator */
BondType& BondType::operator=(const BondType &other)
{
    bond_type = other.bond_type;
    return *this;
}

/** Comparison operator */
bool BondType::operator==(const BondType &other) const
{
    return bond_type == other.bond_type;
}

/** Comparison operator */
bool BondType::operator!=(const BondType &other) const
{
    return not this->operator==(other);
}

const char* BondType::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BondType>() );
}

QString BondType::toString() const
{
    switch (this->bond_type)
    {
        case 0:
            return "undefined";
        case 1:
            return "single";
        case 2:
            return "double";
        case 3:
            return "triple";
        case 4:
            return "aromatic";
        default:
            throw SireError::program_bug(
                QObject::tr("Should not get here: %1").arg(this->bond_type),
                  CODELOC);
    }
}

/** Return the bond type (uses SDF values, e.g. 0 is undefined,
    1 is single, 2 is double, 3 is triple and 4 is aromatic)
*/
int BondType::value() const
{
    return this->bond_type;
}

/** Return the SDF-format value for this bond */
int BondType::sdfValue() const
{
    return this->bond_type;
}

/** Return a single bond */
BondType BondType::singleBond()
{
    return BondType(1);
}

/** Return a double bond */
BondType BondType::doubleBond()
{
    return BondType(2);
}

/** Return a triple bond */
BondType BondType::tripleBond()
{
    return BondType(3);
}

/** Return an aromatic bond */
BondType BondType::aromaticBond()
{
    return BondType(4);
}

/** Return an undefined bond */
BondType BondType::undefinedBond()
{
    return BondType(0);
}

/** Return whether or not the bond type is defined */
bool BondType::isDefined() const
{
    return this->bond_type != 0;
}

/** Return whether or not this is a single bond */
bool BondType::isSingle() const
{
    return this->bond_type == 1;
}

/** Return whether or not this is a double bond */
bool BondType::isDouble() const
{
    return this->bond_type == 2;
}

/** Return whether or not this is a triple bond */
bool BondType::isTriple() const
{
    return this->bond_type == 3;
}

/** Return whether or not this is an aromatic bond */
bool BondType::isAromatic() const
{
    return this->bond_type == 4;
}

