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

#include "radical.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Radical> r_radical;

QDataStream &operator<<(QDataStream &ds, const Radical &r)
{
    writeHeader(ds, r_radical, 1);

    ds << r.radical_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Radical &r)
{
    VersionID v = readHeader(ds, r_radical);

    if (v == 1)
    {
        ds >> r.radical_type;
    }
    else
        throw version_error(v, "1", r_radical, CODELOC);

    return ds;
}

/** Constructor (default is an undefined Radical) */
Radical::Radical()
        : ConcreteProperty<Radical,Property>(), radical_type(0)
{}

/** Construct from the passed string */
Radical::Radical(const QString &str)
        : ConcreteProperty<Radical,Property>()
{
    auto s = str.trimmed().toLower();

    if (s == "singlet")
        this->radical_type = 1;
    else if (s == "doublet")
        this->radical_type = 2;
    else if (s == "triplet")
        this->radical_type = 3;
    else if (s == "undefined")
        this->radical_type = -1;
    else
        throw SireError::invalid_arg(QObject::tr(
            "Cannot interpret radical type '%1'. Should be one of "
            "'singlet', 'doublet', 'triplet' or 'undefined'.")
              .arg(str), CODELOC);
}

/** Construct from the the passed number */
Radical::Radical(int value)
        : ConcreteProperty<Radical,Property>()
{
    if (value >=0 or value <= 3)
    {
        this->radical_type = value;
    }
    else
    {
        throw SireError::invalid_arg(QObject::tr(
            "Invalid radical type '%1'. Should be an integer between "
            "0 to 3").arg(value), CODELOC);
    }
}

/** Copy constructor */
Radical::Radical(const Radical &other)
        : ConcreteProperty<Radical,Property>(other),
          radical_type(other.radical_type)
{}

/** Destructor */
Radical::~Radical()
{}

/** Copy assignment operator */
Radical& Radical::operator=(const Radical &other)
{
    radical_type = other.radical_type;
    return *this;
}

/** Comparison operator */
bool Radical::operator==(const Radical &other) const
{
    return radical_type == other.radical_type;
}

/** Comparison operator */
bool Radical::operator!=(const Radical &other) const
{
    return not this->operator==(other);
}

const char* Radical::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Radical>() );
}

QString Radical::toString() const
{
    switch (this->radical_type)
    {
        case 0:
            return "undefined";
        case 1:
            return "singlet";
        case 2:
            return "doublet";
        case 3:
            return "triplet";
        default:
            throw SireError::program_bug(
                QObject::tr("Should not get here: %1").arg(this->radical_type),
                  CODELOC);
    }
}

/** Return the radical type. 0 is undefined, 1 is singlet, 2 is doublet
    and 3 is triplet
*/
int Radical::value() const
{
    return this->radical_type;
}

/** Return the SDF-format value for this radical. This returns
    '4' if this is a double radical, or 0 otherwise
 */
int Radical::sdfValue() const
{
    if (this->radical_type == 2)
    {
        return 4;
    }
    else
    {
        return 0;
    }
}

/** Return a single Radical */
Radical Radical::singlet()
{
    return Radical(1);
}

/** Return a doublet Radical */
Radical Radical::doublet()
{
    return Radical(2);
}

/** Return a triplet Radical */
Radical Radical::triplet()
{
    return Radical(3);
}

/** Return an undefined Radical */
Radical Radical::undefined()
{
    return Radical(0);
}

/** Return whether or not the Radical is defined */
bool Radical::isDefined() const
{
    return this->radical_type != 0;
}

/** Return whether or not this is a singlet */
bool Radical::isSinglet() const
{
    return this->radical_type == 1;
}

/** Return whether or not this is a doublet */
bool Radical::isDoublet() const
{
    return this->radical_type == 2;
}

/** Return whether or not this is a triplet */
bool Radical::isTriplet() const
{
    return this->radical_type == 3;
}
