/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "pointdipole.h"

#include "SireStream/datastream.h"

using namespace Squire;
using namespace SireStream;

static const RegisterMetaType<PointDipole> r_pointdipole(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const PointDipole &q)
{
    writeHeader(ds, r_pointdipole, 1);
    
    ds << q.cent << q.dipol;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, PointDipole &q)
{
    VersionID v = readHeader(ds, r_pointdipole);
    
    if (v == 1)
    {
        ds >> q.cent >> q.dipol;
    }
    else
        throw version_error(v, "1", r_pointdipole, CODELOC);
        
    return ds;
}

/** Constructor */
PointDipole::PointDipole()
{}

/** Construct a point dipole at the specified location with the 
    specified dipole */
PointDipole::PointDipole(const Vector &coords, const Vector &dipole)
            : cent(coords), dipol(dipole)
{}

/** Copy constructor */
PointDipole::PointDipole(const PointDipole &other)
            : cent(other.cent), dipol(other.dipol)
{}

/** Destructor */
PointDipole::~PointDipole()
{}

/** Copy assignment operator */
PointDipole& PointDipole::operator=(const PointDipole &other)
{
    cent = other.cent;
    dipol = other.dipol;
    return *this;
}

/** Comparison operator */
bool PointDipole::operator==(const PointDipole &other) const
{
    return cent == other.cent and dipol == other.dipol;
}

/** Comparison operator */
bool PointDipole::operator!=(const PointDipole &other) const
{
    return cent != other.cent or dipol != other.dipol;
}

/** Return the location of this point charge */
const Vector& PointDipole::center() const
{
    return cent;
}

/** Return the dipole */
const Vector& PointDipole::dipole() const
{
    return dipol;
}

const char* PointDipole::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PointDipole>() );
}

