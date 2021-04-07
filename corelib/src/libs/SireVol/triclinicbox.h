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

#ifndef SIREVOL_TRICLINICBOX_H
#define SIREVOL_TRICLINICBOX_H

#include "cartesian.h"

#include "SireMaths/matrix.h"

SIRE_BEGIN_HEADER

namespace SireVol
{
class TriclinicBox;
}

SIREVOL_EXPORT QDataStream& operator<<(QDataStream&, const SireVol::TriclinicBox&);
SIREVOL_EXPORT QDataStream& operator>>(QDataStream&, SireVol::TriclinicBox&);

namespace SireVol
{

using SireMaths::Matrix;
using SireMaths::Vector;

/**
A TriclinicBox is a volume  that represents standard periodic boundary conditions
(a 3D box replicated to infinity along all three dimensions).

@author Christopher Woods
*/
class SIREVOL_EXPORT TriclinicBox
        : public SireBase::ConcreteProperty<TriclinicBox,Cartesian>
{

friend SIREVOL_EXPORT QDataStream& ::operator<<(QDataStream&, const TriclinicBox&);
friend SIREVOL_EXPORT QDataStream& ::operator>>(QDataStream&, TriclinicBox&);

public:
    TriclinicBox();
    TriclinicBox(const Vector &v0, const Vector &v1, const Vector &v2);
    TriclinicBox(double a, double b, double c,
                 const SireUnits::Dimension::Angle &alpha,
                 const SireUnits::Dimension::Angle &beta,
                 const SireUnits::Dimension::Angle &gamma);

    TriclinicBox(const TriclinicBox &other);

    ~TriclinicBox();

    TriclinicBox& operator=(const TriclinicBox &other);

    bool operator==(const TriclinicBox &other) const;
    bool operator!=(const TriclinicBox &other) const;

    bool isPeriodic() const;
    bool isCartesian() const;

    QString toString() const;

    /** Get the volume of the triclinic box. */
    SireUnits::Dimension::Volume volume() const;

    /** Set the volume of the triclinic box. */
    SpacePtr setVolume(SireUnits::Dimension::Volume volume) const;

    static const char* typeName();

    /** Whether the triclinic cell has been rotated to comply with the contraints
        of molecular dynamics engines, i.e. vector0 aligned with x axis, vector1
        in x-y plane, and vector2 with positive z component.
     */
    bool isRotated() const;

    double calcDist(const Vector &point0, const Vector &point1) const;

    double calcDist2(const Vector &point0, const Vector &point1) const;

    double calcDist(const CoordGroup &group1, const CoordGroup &group2,
                    DistMatrix &distmat) const;

    double calcDist(const CoordGroup &group, const Vector &point,
                    DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group, const Vector &point,
                     DistMatrix &mat) const;

    double calcDist2(const CoordGroup &group1, const CoordGroup &group2,
                     DistMatrix &distmat) const;

    double calcInvDist(const CoordGroup &group1, const CoordGroup &group2,
                       DistMatrix &distmat) const;

    double calcInvDist2(const CoordGroup &group1, const CoordGroup &group2,
                        DistMatrix &distmat) const;

    DistVector calcDistVector(const Vector &point0, const Vector &point1) const;

    double calcDistVectors(const CoordGroup &group1, const CoordGroup &group2,
                           DistVectorMatrix &distmat) const;

    double calcDistVectors(const CoordGroup &group, const Vector &point,
                           DistVectorMatrix &distmat) const;

    SireUnits::Dimension::Angle calcAngle(const Vector &point0,
                                          const Vector &point1,
                                          const Vector &point2) const;

    SireUnits::Dimension::Angle calcDihedral(const Vector &point0,
                                             const Vector &point1,
                                             const Vector &point2,
                                             const Vector &point3) const;

    bool beyond(double dist, const AABox &aabox0, const AABox &aabox1) const;

    bool beyond(double dist, const CoordGroup &group0,
                const CoordGroup &group1) const;

    double minimumDistance(const CoordGroup &group0, const CoordGroup &group1) const;

    double minimumDistance(const AABox &box0, const AABox &box1) const;

    Vector getMinimumImage(const Vector &point, const Vector &center) const;

    CoordGroup getMinimumImage(const CoordGroup &group, const Vector &center) const;

    CoordGroupArray getMinimumImage(const CoordGroupArray &groups,
                                    const Vector &center,
                                    bool translate_as_one=false) const;

    AABox getMinimumImage(const AABox &aabox, const Vector &center) const;

    QVector<Vector> getImagesWithin(const Vector &point, const Vector &center, double dist) const;

    QList< boost::tuple<double,CoordGroup> >
               getCopiesWithin(const CoordGroup &group,
                               const CoordGroup &center, double dist) const;


    Vector getRandomPoint(const Vector &center, const RanGenerator &generator) const;

	Vector getBoxCenter(const Vector &p) const;
    Vector getBoxCenter(const Vector &p, const Vector &center) const;

    /** Return the first box vector. */
    const Vector& vector0() const;

    /** Return the second box vector. */
    const Vector& vector1() const;

    /** Return the third box vector. */
    const Vector& vector2() const;

    /** Return the rotation matrix. */
    const Matrix& rotationMatrix() const;

    /** Return the cell matrix. */
    Matrix cellMatrix() const;

    /** Return a cubic TriclinicBox with image distance d. */
    static TriclinicBox cubic(double d);

    /** Return a square rhombic dodecahedron TriclinicBox with image distance d. */
    static TriclinicBox rhombicDodecahedronSquare(double d);

    /** Return a hexagonal rhombic dodecahedron TriclinicBox with image distance d. */
    static TriclinicBox rhombicDodecahedronHexagon(double d);

    /** Return a truncated octahedron with image distance d. */
    static TriclinicBox truncatedOctahedron(double d);

protected:

    Vector wrapDelta(const Vector &v0, const Vector &v1) const;

    CoordGroupArray _pvt_getMinimumImage(
                                const CoordGroupArray &groups,
                                const Vector &point) const;

    /** The first box vector */
    Vector v0;

    /** The second box vector */
    Vector v1;

    /** The third box vector */
    Vector v2;

    /** The first box vector (original copy) */
    Vector v0_orig;

    /** The second box vector (original copy) */
    Vector v1_orig;

    /** The third box vector (original copy) */
    Vector v2_orig;

    /** The rotation matrix used to transform the box to meet the requirements
        of molecular dynamics engines.
      */
    Matrix rotation_matrix;

    /** The cell matrix. */
    Matrix cell_matrix;

    /** The inverse of the cell matrix. */
    Matrix cell_matrix_inverse;

    /** The matrix product of cell_matrix_inverse and cell_matrix. */
    Matrix M;

    /** The maximum distance within which a point will always be closer to the
        origin than any of its images.
      */
    double dist_max;

    /** The maximum axis length of the cell. */
    double max_length;

    /** The angle between vectors v0 and v2. */
    double alpha;

    /** The angle between vectors v0 and v2. */
    double beta;

    /** The angle between vectors v1 and v2. */
    double gamma;

    /** The volume of the triclinic cell. */
    double vol;

    /** Whether the triclinic cell has been rotated. */
    bool is_rotated;

    /** The inverse of the lengths of each side of the box */
    Vector invlength;
};

}

Q_DECLARE_METATYPE(SireVol::TriclinicBox)

SIRE_EXPOSE_CLASS( SireVol::TriclinicBox )

SIRE_END_HEADER

#endif
