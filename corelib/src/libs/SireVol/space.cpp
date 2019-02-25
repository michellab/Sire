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

#include "space.h"
#include "cartesian.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

///////////////
/////////////// Implementation of Space
///////////////

static const RegisterMetaType<Space> r_space(MAGIC_ONLY, "SireVol::Space");

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds,
                                       const Space &space)
{
    writeHeader(ds, r_space, 1)
            << static_cast<const Property&>(space);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds,
                                       Space &space)
{
    VersionID v = readHeader(ds, r_space);

    if (v == 1)
    {
        ds >> static_cast<Property&>(space);
    }
    else
        throw version_error(v, "1", r_space, CODELOC);

    return ds;
}

/** Construct a Space. */
Space::Space() : Property()
{}

/** Copy constructor */
Space::Space(const Space &other) : Property(other)
{}

/** Destructor */
Space::~Space()
{}

/** Change the volume of this space by 'delta' */
SpacePtr Space::changeVolume(SireUnits::Dimension::Volume delta) const
{
    return this->setVolume( this->volume() + delta );
}

/** Return a random point within this space using the passed
    random number generator to generate any necessary random
    numbers, and centering the box at the origin */
Vector Space::getRandomPoint(const RanGenerator &generator) const
{
    return this->getRandomPoint(Vector(0,0,0), generator);
}

/** Return a random point within this space using the global 
    random number generator and with the box centered at 'center' */
Vector Space::getRandomPoint(const Vector &center) const
{
    RanGenerator generator;
    return this->getRandomPoint(center, generator);
}

/** Return a random point within this space using the global
    random number generator, and with the box centered at the origin */
Vector Space::getRandomPoint() const
{
    RanGenerator generator;
    return this->getRandomPoint(Vector(0,0,0), generator);
}

/** Assert that 'other' is of the same type as this space

    \throw SireError::incompatible_error
*/
void Space::assertCompatible(const Space &other) const
{
    if ( QLatin1String(this->what()) != QLatin1String(other.what()) )
        throw SireError::incompatible_error( QObject::tr(
            "This space (of type \"%1\") is incompatible with "
            "a space of type \"%2\".")
                .arg(this->what()).arg(other.what()), CODELOC );
}

Q_GLOBAL_STATIC( Cartesian, nullCartesian )

/** Return the default space (Cartesian infinite box) */
const Cartesian& Space::null()
{
    return *(nullCartesian());
}
