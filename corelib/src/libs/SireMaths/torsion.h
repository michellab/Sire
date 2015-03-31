/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMATHS_TORSION_H
#define SIREMATHS_TORSION_H

#include "vector.h"
#include "line.h"
#include "triangle.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Torsion;
}

class QDataStream;
QDataStream& operator<<(QDataStream&, const SireMaths::Torsion&);
QDataStream& operator>>(QDataStream&, SireMaths::Torsion&);

namespace SireMaths
{

using SireUnits::Dimension::Angle;

/**
This class represents a torsion in three dimensional space, e.g. four points 
in space, not necessarily lying in a plane. A torsion is used to calculate 
dihedral angles (imagine each point is an atom). I am not happy with the 
name of this class, and welcome suggestions :-)
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT Torsion
{

friend QDataStream& ::operator<<(QDataStream&, const Torsion&);
friend QDataStream& ::operator>>(QDataStream&, Torsion&);

public:
    Torsion();
    Torsion( const Vector &point0, const Vector &point1,
             const Vector &point2, const Vector &point3 );
    ~Torsion();
    
    static const char* typeName();
    
    const char* what() const
    {
        return Torsion::typeName();
    }

    QString toString() const;

    Angle angle() const;

    Angle improperAngle() const;

    Line line03() const;
    Line line12() const;

    Vector vector03() const;
    Vector vector12() const;

    Triangle triangle1() const;
    Triangle triangle2() const;

    int count() const;

    const Vector& point( int i ) const;
    const Vector& operator[] ( int i ) const;
    const Vector& at( int i ) const;

private:
    /** The four points that make up the torsion */
    Vector points[4];
};

}

Q_DECLARE_METATYPE(SireMaths::Torsion)
Q_DECLARE_TYPEINFO(SireMaths::Torsion, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::Torsion )

SIRE_END_HEADER

#endif
