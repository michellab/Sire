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

#ifndef SIREMATHS_AXISSET_H
#define SIREMATHS_AXISSET_H

#include <QString>
#include <QVector>

#include "matrix.h"
#include "vector.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class AxisSet;
}

class QDataStream;
SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::AxisSet&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::AxisSet&);

namespace SireMaths
{

/**
This class provides a complete set of orthonormal axes that provide a 
frame of reference (origin+axes) for a coordinate system.
 
@author Christopher Woods
*/
class SIREMATHS_EXPORT AxisSet
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const AxisSet&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, AxisSet&);

public:
    AxisSet();
    AxisSet(const Matrix &matrx, Vector orign = Vector());
    AxisSet(const AxisSet &other);
    ~AxisSet();

    static const char* typeName();
    
    const char* what() const
    {
        return AxisSet::typeName();
    }

    QString toString() const;
    
    const Matrix& matrix() const;
    const Matrix& invMatrix() const;
    const Vector& origin() const;

    Vector fromIdentity(const Vector &vec) const;
    QVector<Vector> fromIdentity(const QVector<Vector> &vec) const;
    
    Vector toIdentity(const Vector &vec) const;
    Vector toFrame(const AxisSet &frame, const Vector &vec) const;
    Vector fromFrame(const AxisSet &frame, const Vector &vec) const;
    
    Vector fromIdentity(const Vector &vec, const Vector &delta) const;
    QVector<Vector> fromIdentity(const QVector<Vector> &vecs, const Vector &delta) const;
    
protected:

    /** The matrix that represents the coordinate frame */
    Matrix mat;
    
    /** The inverse of 'mat' */
    Matrix invmat;
    
    /** The origin of this AxisSet */
    Vector orgn;
    
};

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

/** Returns a reference to the matrix representing this AxisSet */
SIRE_ALWAYS_INLINE const Matrix& AxisSet::matrix() const
{
    return mat;
}

/** Returns a reference to the inverse of the matrix representing this AxisSet */
SIRE_ALWAYS_INLINE const Matrix& AxisSet::invMatrix() const
{
    return invmat;
}

/** Returns a reference to the vector representing the origin of this AxisSet */
SIRE_ALWAYS_INLINE const Vector& AxisSet::origin() const
{
    return orgn;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMaths::AxisSet);
Q_DECLARE_TYPEINFO(SireMaths::AxisSet, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMaths::AxisSet )

SIRE_END_HEADER

#endif
