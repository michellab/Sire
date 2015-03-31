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

#include "rational.h"
#include "maths.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<Rational> r_rational(NO_ROOT);

/** Serialise a rational number to a binary datastream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Rational &val)
{
    writeHeader(ds, r_rational, 1) << val.numerator() << val.denominator();
    return ds;
}

/** Deserialise a rational number from a binary datastream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Rational &val)
{
    VersionID v = readHeader(ds, r_rational);

    if (v == 1)
    {
        qint32 num,denom;
        ds >> num >> denom;

        val = Rational(num,denom);
    }
    else
        throw version_error(v, "1", r_rational, CODELOC);

    return ds;
}
