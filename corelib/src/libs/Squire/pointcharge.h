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

#ifndef SQUIRE_POINTCHARGE_H
#define SQUIRE_POINTCHARGE_H

#include "SireMaths/vector.h"
#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class PointCharge;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::PointCharge&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::PointCharge&);

namespace Squire
{

using SireMaths::Vector;

/** This class holds a single point charge. This class is designed
    for speed, and is used within the integral program (the charge
    is held in internal units, and the point is mapped into the
    correct space for the QM program)
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT PointCharge
{

friend QDataStream& ::operator<<(QDataStream&, const PointCharge&);
friend QDataStream& ::operator>>(QDataStream&, PointCharge&);

public:
    PointCharge();
    PointCharge(const Vector &coords, const SireUnits::Dimension::Charge &charge);
    PointCharge(const SireUnits::Dimension::Charge &charge, const Vector &coords);
    
    PointCharge(const PointCharge &other);
    
    ~PointCharge();
    
    static const char* typeName();
    
    PointCharge& operator=(const PointCharge &other);
    
    bool operator==(const PointCharge &other) const;
    bool operator!=(const PointCharge &other) const;
    
    const Vector& center() const;
    double charge() const;
    
private:
    /** The location of this charge, mapped into the correct space */
    Vector cent;
    
    /** The charge, in internal units */
    double q;
};

}

Q_DECLARE_METATYPE( Squire::PointCharge )

SIRE_EXPOSE_CLASS( Squire::PointCharge )

SIRE_END_HEADER

#endif
