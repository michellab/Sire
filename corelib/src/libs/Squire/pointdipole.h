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

#ifndef SQUIRE_POINTDIPOLE_H
#define SQUIRE_POINTDIPOLE_H

#include "SireMaths/vector.h"

SIRE_BEGIN_HEADER

namespace Squire
{
class PointDipole;
}

SQUIRE_EXPORT QDataStream& operator<<(QDataStream&, const Squire::PointDipole&);
SQUIRE_EXPORT QDataStream& operator>>(QDataStream&, Squire::PointDipole&);

namespace Squire
{

using SireMaths::Vector;

/** This class holds a single point dipole. This class is designed
    for speed, and is used within the integral program (the dipole
    is held in internal units, and the point is mapped into the
    correct space for the QM program)
    
    @author Christopher Woods
*/
class SQUIRE_EXPORT PointDipole
{

friend SQUIRE_EXPORT QDataStream& ::operator<<(QDataStream&, const PointDipole&);
friend SQUIRE_EXPORT QDataStream& ::operator>>(QDataStream&, PointDipole&);

public:
    PointDipole();
    PointDipole(const Vector &coords, const Vector &dipole);
    
    PointDipole(const PointDipole &other);
    
    ~PointDipole();
    
    static const char* typeName();
    
    PointDipole& operator=(const PointDipole &other);
    
    bool operator==(const PointDipole &other) const;
    bool operator!=(const PointDipole &other) const;
    
    const Vector& center() const;
    const Vector& dipole() const;
    
private:
    /** The location of this dipole, mapped into the correct space */
    Vector cent;
    
    /** The dipole, in internal units */
    Vector dipol;
};

}

Q_DECLARE_METATYPE( Squire::PointDipole )

SIRE_EXPOSE_CLASS( Squire::PointDipole )

SIRE_END_HEADER

#endif
