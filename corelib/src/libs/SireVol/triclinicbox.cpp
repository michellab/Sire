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

#include <QDebug>
#include <QList>
#include <QMutex>
#include <QPair>

#include <limits>
#include <cmath>

#ifdef SIRE_USE_SSE
    #ifdef __SSE__
        #include <emmintrin.h>   // CONDITIONAL_INCLUDE
    #else
        #undef SIRE_USE_SSE
    #endif
#endif

#include "triclinicbox.h"
#include "coordgroup.h"

#include "SireMaths/rangenerator.h"

#include "SireBase/countflops.h"

#include "SireError/errors.h"
#include "SireStream/datastream.h"

#include <QDebug>

// Helper struct for sorting based on first pair value.
struct QPairFirstComparer
{
    template<typename T1, typename T2>
    bool operator()(const QPair<T1,T2> & a, const QPair<T1,T2> & b) const
    {
        return a.first < b.first;
    }
};

using namespace SireVol;
using namespace SireBase;
using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<TriclinicBox> r_box;

/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const TriclinicBox &box)
{
    writeHeader(ds, r_box, 1)
               << box.v0_orig
               << box.v1_orig
               << box.v2_orig;

               //no need to store anything else as it can be regenerated

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream &operator>>(QDataStream &ds, TriclinicBox &box)
{
    VersionID v = readHeader(ds, r_box);

    if (v == 1)
    {
        Vector v0, v1, v2;

        // Load in the original box vectors.
        ds >> v0 >> v1 >> v2;

        // Reconstruct the box.
        box = TriclinicBox(v0, v1, v2);
    }
    else
        throw version_error(v, "1", r_box, CODELOC);

    return ds;
}

/** This is the maximum dimension of the box (so that .volume() doesn't overflow) */
const double max_boxlength = std::pow(0.9 * std::numeric_limits<double>::max(), 1.0/3.0);

/** Construct a default TriclinicBox volume (maximum volume) */
TriclinicBox::TriclinicBox() : ConcreteProperty<TriclinicBox,Cartesian>()
{
    Vector v0(max_boxlength, 0, 0);
    Vector v1(0, max_boxlength, 0);
    Vector v2(0, 0, max_boxlength);

    *this = TriclinicBox(v0, v1, v2);
}

/** Construct a TriclinicBox with the specified lattice vectors */
TriclinicBox::TriclinicBox(const Vector &v0,
                           const Vector &v1,
                           const Vector &v2)
            : ConcreteProperty<TriclinicBox,Cartesian>(),
              v0(v0),
              v1(v1),
              v2(v2),
              v0_orig(v0),
              v1_orig(v1),
              v2_orig(v2)
{
    // What follows was adapted from boxmtx.pl by Tsjerk A. Wassenaar.

    // Get the magnitudes of the lattice box vectors.
    auto m0 = v0.magnitude();
    auto m1 = v1.magnitude();
    auto m2 = v2.magnitude();

    // Create a list containing of magnitude/vector pairs.
    QList<QPair<double, Vector> > magnitudes;
    magnitudes.append(qMakePair(m0, v0));
    magnitudes.append(qMakePair(m1, v1));
    magnitudes.append(qMakePair(m2, v2));

    // Sort the pairs in ascending order.
    qSort(magnitudes.begin(), magnitudes.end(), QPairFirstComparer());

    // Reorder the vectors in order of descending magnitude.
    this->v0 = magnitudes[2].second;
    this->v1 = magnitudes[1].second;
    this->v2 = magnitudes[0].second;

    // Project v1 onto v0.
    auto xn1 = Vector::dot(v0, v1) / m0;

    // Evaluate the square root of the difference between the length of v1 squared
    // and the projection of that vector on v0 squared.
    auto yn1 = sqrt(m1*m1 - xn1*xn1);

    // Evaluate the vector product of v0 and v1 and take the magnitude.
    auto zn2 = Vector::realCross(v0, v1).magnitude();

    // Construct the original matrix.
    Matrix X = Matrix(v0, v1, Vector::realCross(v0, v1)).transpose();

    // Constuct the rotated matrix.
    Matrix X_rot(Vector(m0, xn1, 0),
                 Vector(0, yn1, 0),
                 Vector(0, 0, zn2));

    // Evaluate and store the rotation matrix and its inverse.
    this->rotation_matrix = X_rot * X.inverse();

    // Rotate the lattice vectors.
    this->v0 = this->rotation_matrix * v0;
    this->v1 = this->rotation_matrix * v1;
    this->v2 = this->rotation_matrix * v2;

    // Make sure v2 has a positive z component.
    if (this->v2.z() < 0)
        this->v2 *= -1;

    /* Next perform a lattice reduction such that the following conditions are
       met:

       | v1.x | <= 1/2 v0.x
       | v2.x | <= 1/2 v0.x
       | v2.y | <= 1/2 v1.y

       These constraints can be solved using simultaneous equations.
     */

    // v1.x
    double v1x;
    if (this->v1.x() != 0)
    {
        v1x = this->v1.x()
            - int((this->v1.x() / this->v0.x())
            + (this->v1.x() / (2.0 * qAbs(this->v1.x()))))
            * this->v0.x();
    }
    else
    {
        v1x = this->v1.x();
    }

    // v2.x
    double v2x;
    if (this->v2.x() != 0)
    {
        v2x = this->v2.x()
            - int((this->v2.x() / this->v1.x())
            + (this->v2.x() / (2.0 * qAbs(this->v2.x()))))
            * this->v1.x();
    }
    else
    {
        v2x = this->v2.x();
    }

    // v2.y
    double v2y;
    if (this->v2.y() != 0)
    {
        v2y = this->v2.y()
            - int((this->v2.y() / this->v1.y())
            + (this->v2.y() / (2.0 * qAbs(this->v2.y()))))
            * this->v1.y();
    }
    else
    {
        v2y = this->v2.y();
    }

    // Update the lattice vectors.
    this->v1 = Vector(v1x, this->v1.y(), this->v1.z());
    this->v2 = Vector(v2x, v2y, this->v2.z());

    // Store the cell matrix and its inverse.
    this->cell_matrix = Matrix(this->v0, this->v1, this->v2).transpose();
    this->cell_matrix_inverse = this->cell_matrix.inverse();

    // Store the product of the transpose of cell_matrix and cell_matrix.
    this->M = this->cell_matrix.transpose() * this->cell_matrix;

    // Work out the maximum distance for minimum image calculations.

    // Get the magnitudes of the updated lattice box vectors. (Half-lengths)
    m0 = this->v0.magnitude() * 0.5;
    m1 = this->v1.magnitude() * 0.5;
    m2 = this->v2.magnitude() * 0.5;

    // Store the minimum half distance.
    if ((m0 < m1) and (m0 < m2))
    {
        this->dist_max = m0;
    }
    else if ((m0 < m1) and (m2 < m0))
    {
        this->dist_max = m2;
    }
    else if (m1 < m2)
    {
        this->dist_max = m1;
    }
    else
    {
        this->dist_max = m2;
    }

    // Scale up half-lengths again.
    m0 *= 2.0;
    m1 *= 2.0;
    m2 *= 2.0;

    // Work out Ghe angle between each pair of vectors.
    this->alpha = Vector::angle(this->v1, this->v2).value();
    this->beta = Vector::angle(this->v0, this->v2).value();
    this->gamma = Vector::angle(this->v1, this->v0).value();

    auto cos_alpha = cos(this->alpha);
    auto cos_beta = cos(this->beta);
    auto cos_gamma = cos(this->gamma);

    // Now work out the volume of the cell.
    this->vol = this->v0.magnitude() * this->v1.magnitude() * this->v2.magnitude() *
                std::sqrt(1 - cos_alpha*cos_alpha
                            - cos_beta*cos_beta
                            - cos_gamma*cos_gamma
                            + 2.0*cos_alpha*cos_beta*cos_gamma);
}

/** Copy constructor */
TriclinicBox::TriclinicBox(const TriclinicBox &other)
            : ConcreteProperty<TriclinicBox,Cartesian>(other),
              v0(other.v0),
              v1(other.v1),
              v2(other.v2),
              v0_orig(other.v0_orig),
              v1_orig(other.v1_orig),
              v2_orig(other.v2_orig),
              rotation_matrix(other.rotation_matrix),
              cell_matrix(other.cell_matrix),
              cell_matrix_inverse(other.cell_matrix_inverse),
              M(other.M),
              dist_max(other.dist_max),
              alpha(other.alpha),
              beta(other.beta),
              gamma(other.gamma),
              vol(other.vol)
{}

/** Destructor */
TriclinicBox::~TriclinicBox()
{}

/** Copy assignment operator */
TriclinicBox& TriclinicBox::operator=(const TriclinicBox &other)
{
    if (this != &other)
    {
        v0 = other.v0;
        v1 = other.v1;
        v2 = other.v2;
        v0_orig = other.v0_orig;
        v1_orig = other.v1_orig;
        v2_orig = other.v2_orig;
        rotation_matrix = other.rotation_matrix;
        cell_matrix = other.cell_matrix;
        dist_max = other.dist_max;
        alpha = other.alpha;
        beta = other.beta;
        gamma = other.gamma;
        vol = other.vol;
        Cartesian::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool TriclinicBox::operator==(const TriclinicBox &other) const
{
    return v0 == other.v0 and
           v1 == other.v1 and
           v2 == other.v2 and
           v0_orig == other.v0_orig and
           v1_orig == other.v1_orig and
           v2_orig == other.v2_orig;
}

/** Comparison operator */
bool TriclinicBox::operator!=(const TriclinicBox &other) const
{
    return v0 != other.v0 or
           v1 != other.v1 or
           v2 != other.v2 or
           v0_orig != other.v0_orig or
           v1_orig != other.v1_orig or
           v2_orig != other.v2_orig;
}

/** A Triclinic box is periodic! */
bool TriclinicBox::isPeriodic() const
{
    return true;
}

/** A Triclinic box is cartesian */
bool TriclinicBox::isCartesian() const
{
    return true;
}

/** Return the first box vector */
const Vector& TriclinicBox::vector0() const
{
    return v0;
}

/** Return the second box vector */
const Vector& TriclinicBox::vector1() const
{
    return v1;
}

/** Return the third box vector */
const Vector& TriclinicBox::vector2() const
{
    return v2;
}

/** Return the rotation matrix */
const Matrix& TriclinicBox::rotationMatrix() const
{
    return rotation_matrix;
}

/** Return the cell matrix */
Matrix TriclinicBox::cellMatrix() const
{
    Matrix M = Matrix(v0, v1, v2).transpose();
    return M;
}

/** Return a cubic TriclinicBox with image distance d. */
TriclinicBox TriclinicBox::cubic(double d)
{
    Vector v0(d, 0, 0);
    Vector v1(0, d, 0);
    Vector v2(0, 0, d);

    return TriclinicBox(v0, v1, v2);
}

/** Return a square rhombic dodecahedron TriclinicBox with image distance d. */
TriclinicBox TriclinicBox::rhombicDodecahedronSquare(double d)
{
    Vector v0(d, 0, 0);
    Vector v1(0, d, 0);
    Vector v2(0.5*d, 0.5*d, 0.5*sqrt(2)*d);

    return TriclinicBox(v0, v1, v2);
}

/** Return a hexagonal rhombic dodecahedron TriclinicBox with image distance d. */
TriclinicBox TriclinicBox::rhombicDodecahedronHexagon(double d)
{
    Vector v0(d, 0, 0);
    Vector v1(0.5, 0.5*sqrt(3)*d, 0);
    Vector v2(0.5*d, (1/6.0)*sqrt(3)*d, (1/3.0)*sqrt(6)*d);

    return TriclinicBox(v0, v1, v2);
}

/** Return a truncated octahedron with image distance d. */
TriclinicBox TriclinicBox::truncatedOctahedron(double d)
{
    Vector v0(d, 0, 0);
    Vector v1(d/3.0, (2/3.0)*sqrt(2)*d, 0);
    Vector v2(-d/3.0, (1/3.0)*sqrt(2)*d, (1/3.0)*sqrt(6)*d);

    return TriclinicBox(v0, v1, v2);
}

/** Return a string representation of this space */
QString TriclinicBox::toString() const
{
    return QObject::tr("TriclinicBox( %1, %2, %3 )")
                .arg( this->vector0().toString() )
                .arg( this->vector1().toString() )
                .arg( this->vector2().toString() );
}

/** Return the volume of the central box of this space. */
SireUnits::Dimension::Volume TriclinicBox::volume() const
{
    return SireUnits::Dimension::Volume(this->vol);
}

/** Calculate the delta that needs to be subtracted from the interatomic
    distances so that the molecules are all wrapped into the same triclinic box */
Vector TriclinicBox::wrapDelta(const Vector &v0, const Vector &v1) const
{
    // Work out the distance vector in "box" space.
    auto dist_box = this->cell_matrix_inverse*v1 - this->cell_matrix_inverse*v0;

    // Return the shifts over the three lattice vectors.
    return Vector(int(dist_box.x()), int(dist_box.y()), int(dist_box.z()));
}

/** Return a copy of this space with the volume of set to 'volume'
    - this will scale the space uniformly, keeping the center at
    the same location, to achieve this volume */
SpacePtr TriclinicBox::setVolume(SireUnits::Dimension::Volume vol) const
{
    double old_volume = this->volume();
    double new_volume = vol;

    if (new_volume < 0)
        throw SireError::invalid_arg( QObject::tr(
            "You cannot set the volume of a periodic box to a negative value! (%1)")
                .arg(new_volume), CODELOC );

    if (old_volume == new_volume)
        return *this;

    double scl = std::pow( new_volume / old_volume, 1.0/3.0 ); // rats - no cbrt function!

    return TriclinicBox(scl*this->v0, scl*this->v1, scl*this->v2);
}

/** Calculate the distance between two points */
double TriclinicBox::calcDist(const Vector &point0, const Vector &point1) const
{
    // Work out the distance vector in "box" space.
    auto dist_box = this->cell_matrix_inverse*point1 - this->cell_matrix_inverse*point0;

    // Extract the fractional components of the distance.
    double x = dist_box.x() - int(dist_box.x());
    double y = dist_box.y() - int(dist_box.y());
    double z = dist_box.z() - int(dist_box.z());

    // Shift to box.
    if (x >= 0.5) x -= 1.0;
    if (y >= 0.5) y -= 1.0;
    if (z >= 0.5) z -= 1.0;

    // Construct a vector from the fractional components.
    Vector fract_dist(x, y, z);

    // Calculate the distance in the minimum image.
    auto dist = std::sqrt(Vector::dot(fract_dist, this->M*fract_dist));

    return dist;
}

/** Calculate the distance squared between two points */
double TriclinicBox::calcDist2(const Vector &point0, const Vector &point1) const
{
    auto dist = this->calcDist(point0, point1);
    return dist*dist;
}

/** Populate the matrix 'mat' with the distances between all of the
    atoms of the two CoordGroups. Return the shortest distance^2 between the two
    CoordGroups. */
double TriclinicBox::calcDist(const CoordGroup &group0, const CoordGroup &group1,
                              DistMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            const double dist = this->calcDist(point0, array1[j]);

            mindist = qMin(mindist, dist);
            mat[j] = dist;
        }
    }

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'mat' with the distances between all of the
    atoms of the passed CoordGroup to the passed point. Return the shortest
    distance. */
double TriclinicBox::calcDist(const CoordGroup &group, const Vector &point,
                              DistMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointer to the array - this provides more efficient access
    const Vector *array = group.constData();

    mat.setOuterIndex(0);

    for (int j=0; j<n; ++j)
    {
        const double dist = this->calcDist(point, array[j]);

        mindist = qMin(mindist, dist);
        mat[j] = dist;
    }

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'mat' with the distances squared between all of the
    atoms of the passed CoordGroup to the passed point. Return the shortest
    distance. */
double TriclinicBox::calcDist2(const CoordGroup &group, const Vector &point,
                               DistMatrix &mat) const
{
    double mindist2(std::numeric_limits<double>::max());

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointer to the array - this provides more efficient access
    const Vector *array = group.constData();

    mat.setOuterIndex(0);

    for (int j=0; j<n; ++j)
    {
        const double dist2 = this->calcDist2(point, array[j]);

        mindist2 = qMin(mindist2, dist2);
        mat[j] = dist2;
    }

    //return the minimum distance
    return sqrt(mindist2);
}

/** Populate the matrix 'mat' with the distances^2 between all of the
    atoms of the two CoordGroups. Return the shortest distance between the
    two CoordGroups. */
double TriclinicBox::calcDist2(const CoordGroup &group0, const CoordGroup &group1,
                               DistMatrix &mat) const
{
    double mindist2(std::numeric_limits<double>::max());

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two atoms
            const double tmpdist = this->calcDist(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist2 = qMin(tmpdist,mindist2);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    //return the minimum distance
    return sqrt(mindist2);
}

/** Populate the matrix 'mat' with the inverse distances between all of the
    atoms of the two CoordGroups. Return the shortest distance between the two CoordGroups. */
double TriclinicBox::calcInvDist(const CoordGroup &group0, const CoordGroup &group1,
                                 DistMatrix &mat) const
{
    double maxinvdist(0);
    double tmpdist;

    int n0 = group0.count();
    int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two atoms
            tmpdist = 1.0 / this->calcDist(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist = qMax(tmpdist,maxinvdist);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    //return the shortest distance
    return 1.0 / maxinvdist;
}

/** Populate the matrix 'mat' with the inverse distances^2 between all of the
    atoms of the two CoordGroups. Return the shortest distance between the two CoordGroups. */
double TriclinicBox::calcInvDist2(const CoordGroup &group0, const CoordGroup &group1,
                                  DistMatrix &mat) const
{
    double maxinvdist2(0);
    double tmpdist;

    int n0 = group0.count();
    int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two atoms
            tmpdist = 1.0 / this->calcDist(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist2 = qMax(tmpdist,maxinvdist2);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    //return the shortest distance
    return 1.0 / sqrt(maxinvdist2);
}

/** Calculate the distance vector between two points */
DistVector TriclinicBox::calcDistVector(const Vector &point0,
                                        const Vector &point1) const
{
    // Work out the distance vector in "box" space.
    auto dist_box = this->cell_matrix_inverse*point0 - this->cell_matrix_inverse*point1;

    // Extract the fractional components of the distance.
    double x = dist_box.x() - int(dist_box.x());
    double y = dist_box.y() - int(dist_box.y());
    double z = dist_box.z() - int(dist_box.z());

    // Shift to box.
    if (x >= 0.5) x -= 1.0;
    if (y >= 0.5) y -= 1.0;
    if (z >= 0.5) z -= 1.0;

    // Construct a vector from the fractional components.
    Vector fract_dist(x, y, z);

    // Return the fractional distance vector mapped back to the triclinic
    // cell space.
    return this->cellMatrix()*fract_dist;
}

/** Populate the matrix 'distmat' between all the points of the two CoordGroups
    'group1' and 'group2' - the returned matrix has the vectors pointing
    from each point in 'group1' to each point in 'group2'. This returns
    the shortest distance between two points in the group */
double TriclinicBox::calcDistVectors(const CoordGroup &group0, const CoordGroup &group1,
                                     DistVectorMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        //add the delta to the coordinates of atom0
        Vector point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            mat[j] = this->calcDistVector(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist = qMin(mat[j].length(),mindist);
        }
    }

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'distmat' between all the points passed CoordGroup
    to the point 'point' - the returned matrix has the vectors pointing
    from the point to each point in 'group'. This returns
    the shortest distance. */
double TriclinicBox::calcDistVectors(const CoordGroup &group, const Vector &point,
                                     DistVectorMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array = group.constData();

    mat.setOuterIndex(0);

    for (int j=0; j<n; ++j)
    {
        mat[j] = this->calcDistVector(point, array[j]);

        //store the minimum distance, the value expected to be the minimum
        //value is most efficiently placed as the second argument
        mindist = qMin(mat[j].length(),mindist);
    }

    //return the minimum distance
    return mindist;
}

const char* TriclinicBox::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TriclinicBox>() );
}