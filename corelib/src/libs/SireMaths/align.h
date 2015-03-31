/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2014  Christopher Woods
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

#ifndef SIREMATHS_ALIGN_H
#define SIREMATHS_ALIGN_H

#include "SireMaths/vector.h"
#include "SireMaths/matrix.h"
#include "SireMaths/quaternion.h"
#include "SireMaths/axisset.h"

#include <QVector>
#include <QList>

SIRE_BEGIN_HEADER

namespace SireMaths
{
class Transform;
}

QDataStream& operator<<(QDataStream&, const SireMaths::Transform&);
QDataStream& operator>>(QDataStream&, SireMaths::Transform&);

namespace SireMaths
{
    Vector getCentroid(const QVector<Vector> &p, int n=-1);
    
    double getRMSD(const QVector<Vector> &p, const QVector<Vector> &q, int n=-1);

    /** This class holds everything needed to apply a transformation on a set
        of points. This holds the amount by which to translate the points, together
        with the center of rotation and amount by which to rotate
        
        @author Christopher Woods
    */
    class SIREMATHS_EXPORT Transform
    {
    
    friend QDataStream& ::operator<<(QDataStream&, const Transform&);
    friend QDataStream& ::operator>>(QDataStream&, Transform&);
    
    public:
        Transform();
        Transform(const Vector &delta);
        
        Transform(const Quaternion &rotmat, const Vector &center);
        Transform(const Matrix &rotmat, const Vector &center);
        
        Transform(const Vector &delta, const Quaternion &rotmat, const Vector &center);
        Transform(const Vector &delta, const Matrix &rotmat, const Vector &center);
        
        Transform(const Transform &other);
        
        ~Transform();
        
        Transform& operator=(const Transform &other);
        
        bool operator==(const Transform &other) const;
        bool operator!=(const Transform &other) const;
        
        static const char* typeName();
        const char* what() const;
        
        QString toString() const;
        
        bool isNull() const;
        bool isZero() const;
        
        Vector operator()(const Vector &point) const;
        QVector<Vector> operator()(const QVector<Vector> &point) const;
        
        Vector apply(const Vector &point) const;
        QVector<Vector> apply(const QVector<Vector> &points) const;
        
        Vector* apply(Vector *coords, int sz) const;
        
        Vector translationDelta() const;
        
        Vector rotationCenter() const;
        
        Quaternion rotationQuaternion() const;
        Matrix rotationMatrix() const;
        
    private:
        /** Amount by which to translate */
        Vector delta;
        
        /** The center of rotation */
        Vector rotcent;
        
        /** The rotation matrix (as a quaternion) */
        Quaternion rotmat;
    };

    Transform kabaschFit(const QVector<Vector> &p,
                         const QVector<Vector> &q);

    Matrix kabasch(const QVector<Vector> &p,
                   const QVector<Vector> &q);

    Transform getAlignment(const QVector<Vector> &p,
                           const QVector<Vector> &q,
                           bool fit=true);
    
    QVector<Vector> align(const QVector<Vector> &p,
                          const QVector<Vector> &q,
                          bool fit=true);
}

Q_DECLARE_METATYPE( SireMaths::Transform )

SIRE_EXPOSE_CLASS( SireMaths::Transform )

SIRE_EXPOSE_FUNCTION( SireMaths::getCentroid )
SIRE_EXPOSE_FUNCTION( SireMaths::getRMSD )
SIRE_EXPOSE_FUNCTION( SireMaths::kabasch )
SIRE_EXPOSE_FUNCTION( SireMaths::kabaschFit )
SIRE_EXPOSE_FUNCTION( SireMaths::getAlignment )
SIRE_EXPOSE_FUNCTION( SireMaths::align )

SIRE_END_HEADER

#endif
