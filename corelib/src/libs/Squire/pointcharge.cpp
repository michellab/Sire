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

#include "pointcharge.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace Squire;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

static const RegisterMetaType<PointCharge> r_pointcharge(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const PointCharge &q)
{
    writeHeader(ds, r_pointcharge, 1);
    
    ds << q.cent << q.q;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, PointCharge &q)
{
    VersionID v = readHeader(ds, r_pointcharge);
    
    if (v == 1)
    {
        ds >> q.cent >> q.q;
    }
    else
        throw version_error(v, "1", r_pointcharge, CODELOC);
        
    return ds;
}

/** Constructor */
PointCharge::PointCharge() : q(0)
{}

/** Construct a point charge at the specified location with the 
    specified charge */
PointCharge::PointCharge(const Vector &coords, const Charge &charge)
            : cent(coords), q(charge.to(mod_electron))
{}

/** Construct a point charge at the specified location with the 
    specified charge */
PointCharge::PointCharge(const Charge &charge, const Vector &coords)
            : cent(coords), q(charge.to(mod_electron))
{}

/** Copy constructor */
PointCharge::PointCharge(const PointCharge &other)
            : cent(other.cent), q(other.q)
{}

/** Destructor */
PointCharge::~PointCharge()
{}

/** Copy assignment operator */
PointCharge& PointCharge::operator=(const PointCharge &other)
{
    cent = other.cent;
    q = other.q;
    return *this;
}

/** Comparison operator */
bool PointCharge::operator==(const PointCharge &other) const
{
    return q == other.q and cent == other.cent;
}

/** Comparison operator */
bool PointCharge::operator!=(const PointCharge &other) const
{
    return q != other.q or cent != other.cent;
}

/** Return the location of this point charge */
const Vector& PointCharge::center() const
{
    return cent;
}

/** Return the magnitude of this point charge (in internal units) */
double PointCharge::charge() const
{
    return q;
}

const char* PointCharge::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PointCharge>() );
}
