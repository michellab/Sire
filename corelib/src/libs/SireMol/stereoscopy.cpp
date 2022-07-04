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

#include "stereoscopy.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<Stereoscopy> r_stereo;

QDataStream &operator<<(QDataStream &ds, const Stereoscopy &s)
{
    writeHeader(ds, r_stereo, 1);

    ds << s.stereo_type;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, Stereoscopy &s)
{
    VersionID v = readHeader(ds, r_stereo);

    if (v == 1)
    {
        ds >> s.stereo_type;
    }
    else
        throw version_error(v, "1", r_stereo, CODELOC);

    return ds;
}

/** Constructor (default is an undefined stereoscopy) */
Stereoscopy::Stereoscopy()
            : ConcreteProperty<Stereoscopy,Property>(), stereo_type(-1)
{}

/** Construct from the passed string */
Stereoscopy::Stereoscopy(const QString &str)
            : ConcreteProperty<Stereoscopy,Property>()
{
    auto s = str.trimmed().toLower();

    if (s == "up")
        this->stereo_type = 1;
    else if (s == "down")
        this->stereo_type = 6;
    else if (s == "not stereo")
        this->stereo_type = 0;
    else if (s == "undefined")
        this->stereo_type = -1;
    else
        throw SireError::invalid_arg(QObject::tr(
            "Cannot interpret stereo type '%1'. Should be one of "
            "'up', 'down', 'not stereo' or 'undefined'.")
              .arg(str), CODELOC);
}

/** Construct from the the passed number */
Stereoscopy::Stereoscopy(int value)
            : ConcreteProperty<Stereoscopy,Property>()
{
    if (value == 0 or value == 1 or value == -1 or value == 6)
    {
        this->stereo_type = value;
    }
    else
    {
        throw SireError::invalid_arg(QObject::tr(
            "Invalid stereo type '%1'. Should be an integer in "
            "[-1, 0, 1, 6]").arg(value), CODELOC);
    }
}

/** Copy constructor */
Stereoscopy::Stereoscopy(const Stereoscopy &other)
         : ConcreteProperty<Stereoscopy,Property>(other),
           stereo_type(other.stereo_type)
{}

/** Destructor */
Stereoscopy::~Stereoscopy()
{}

/** Copy assignment operator */
Stereoscopy& Stereoscopy::operator=(const Stereoscopy &other)
{
    stereo_type = other.stereo_type;
    return *this;
}

/** Comparison operator */
bool Stereoscopy::operator==(const Stereoscopy &other) const
{
    return stereo_type == other.stereo_type;
}

/** Comparison operator */
bool Stereoscopy::operator!=(const Stereoscopy &other) const
{
    return not this->operator==(other);
}

const char* Stereoscopy::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Stereoscopy>() );
}

QString Stereoscopy::toString() const
{
    switch (this->stereo_type)
    {
        case 0:
            return "not stereo";
        case 1:
            return "up";
        case 6:
            return "down";
        case -1:
            return "undefined";
        default:
            throw SireError::program_bug(
                QObject::tr("Should not get here: %1").arg(this->stereo_type),
                  CODELOC);
    }
}

/** Return the stereo type (uses SDF values, e.g. 0 is not stereo,
    1 is up, 6 is down. We have added -1 to mean undefined)
*/
int Stereoscopy::value() const
{
    return this->stereo_type;
}

/** Return the SDF-format value for this bond. This returns
    '0' if the stereoscopy is undefined
 */
int Stereoscopy::sdfValue() const
{
    if (this->stereo_type == -1)
    {
        return 0;
    }
    else
    {
        return this->stereo_type;
    }
}

/** Return an up Stereoscopy */
Stereoscopy Stereoscopy::up()
{
    return Stereoscopy(1);
}

/** Return a down Stereoscopy */
Stereoscopy Stereoscopy::down()
{
    return Stereoscopy(6);
}

/** Return a "not stereo" Stereoscopy */
Stereoscopy Stereoscopy::notStereo()
{
    return Stereoscopy(0);
}

/** Return an undefined Stereoscopy */
Stereoscopy Stereoscopy::undefined()
{
    return Stereoscopy(-1);
}

/** Return whether or not the stereoscopy is defined */
bool Stereoscopy::isDefined() const
{
    return this->stereo_type != -1;
}

/** Return whether or not this is an up bond */
bool Stereoscopy::isUp() const
{
    return this->stereo_type == 1;
}

/** Return whether or not this is a down bond */
bool Stereoscopy::isDown() const
{
    return this->stereo_type == 6;
}

/** Return whether or not this is a "not stereo" bond */
bool Stereoscopy::isNotStereo() const
{
    return this->stereo_type == 0;
}
