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
            - int((this->v2.x() / this->v0.x())
            + (this->v2.x() / (2.0 * qAbs(this->v2.x()))))
            * this->v0.x();
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

    // Store the maximum length.
    if ((m0 > m1) and (m0 > m2))
    {
        this->max_length = m0;
    }
    else if ((m0 > m1) and (m2 > m0))
    {
        this->max_length = m2;
    }
    else if (m1 > m2)
    {
        this->max_length = m1;
    }
    else
    {
        this->max_length = m2;
    }

    // Work out the angle between each pair of vectors.
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

    // Store the inverse lengths of the lattice vectors.
    this->invlength = Vector(1.0/this->v0.magnitude(),
                             1.0/this->v1.magnitude(),
                             1.0/this->v2.magnitude());
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
              max_length(other.max_length),
              alpha(other.alpha),
              beta(other.beta),
              gamma(other.gamma),
              vol(other.vol),
              invlength(other.invlength)
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
        max_length = other.max_length;
        alpha = other.alpha;
        beta = other.beta;
        gamma = other.gamma;
        vol = other.vol;
        invlength = other.invlength;
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

/** In general, a triclinic box isn't Cartesian. **/
bool TriclinicBox::isCartesian() const
{
    // Only cubic boxes are Cartesian.
    return this->alpha == M_PI_2 and
           this->beta  == M_PI_2 and
           this->gamma == M_PI_2;
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

/** Calculate the delta that needs to be subtracted from the interatomic
    distances so that the molecules are all wrapped into the same periodic box */
Vector TriclinicBox::wrapDelta(const Vector &v0, const Vector &v1) const
{
    // Work out the positions of v0 and v1 in "box" space.
    auto v0_box = this->cell_matrix_inverse*v0;
    auto v1_box = this->cell_matrix_inverse*v1;

    // Work out the distance vector in "box" space.
    auto dist_box = v1_box - v0_box;

    // Extract the integer and fractional components of the distance.
    int int_x = int(dist_box.x());
    int int_y = int(dist_box.y());
    int int_z = int(dist_box.z());
    double frac_x = dist_box.x() - int_x;
    double frac_y = dist_box.y() - int_y;
    double frac_z = dist_box.z() - int_z;

    // Shift to box.

    // x
    if      (frac_x >=  0.5) int_x += 1.0;
    else if (frac_x <= -0.5) int_x -= 1.0;

    // y
    if      (frac_y >=  0.5) int_y += 1.0;
    else if (frac_y <= -0.5) int_y -= 1.0;

    // z
    if      (frac_z >=  0.5) int_z += 1.0;
    else if (frac_z <= -0.5) int_z -= 1.0;

    // Return the shifts over the box vectors.
    return this->cell_matrix*Vector(int_x, int_y, int_z);
}

/** Calculate the distance between two points */
double TriclinicBox::calcDist(const Vector &point0, const Vector &point1) const
{
    // Work out the distance vector in "box" space.
    auto dist_box = this->cell_matrix_inverse*point1 - this->cell_matrix_inverse*point0;

    // Extract the integer and fractional components of the distance.
    int int_x = int(dist_box.x());
    int int_y = int(dist_box.y());
    int int_z = int(dist_box.z());
    double frac_x = dist_box.x() - int_x;
    double frac_y = dist_box.y() - int_y;
    double frac_z = dist_box.z() - int_z;

    // Shift to box.

    // x
    if      (frac_x >=  0.5) frac_x -= 1.0;
    else if (frac_x <= -0.5) frac_x += 1.0;

    // y
    if      (frac_y >=  0.5) frac_y -= 1.0;
    else if (frac_y <= -0.5) frac_y += 1.0;

    // z
    if      (frac_z >=  0.5) frac_z -= 1.0;
    else if (frac_z <= -0.5) frac_z += 1.0;

    // Construct a vector from the fractional components.
    Vector frac_dist(frac_x, frac_y, frac_z);

    // Calculate the distance in the minimum image.
    auto dist = std::sqrt(Vector::dot(frac_dist, this->M*frac_dist));

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i] + wrapdelta;
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            const double dist = Vector::distance(point0, array1[j]);

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group.aaBox().center(), point);
    Vector wrapped_point = point + wrapdelta;

    for (int j=0; j<n; ++j)
    {
        const double dist = Vector::distance(wrapped_point, array[j]);

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group.aaBox().center(), point);
    Vector wrapped_point = point + wrapdelta;

    mat.setOuterIndex(0);

    for (int j=0; j<n; ++j)
    {
        const double dist2 = Vector::distance2(wrapped_point, array[j]);

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i] + wrapdelta;
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            const double tmpdist = Vector::distance2(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist2 = qMin(tmpdist, mindist2);

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i] + wrapdelta;
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            tmpdist = Vector::invDistance(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist = qMax(tmpdist, maxinvdist);

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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i] + wrapdelta;
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the inverse squared distance between the two points
            tmpdist = Vector::invDistance2(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist2 = qMax(tmpdist, maxinvdist2);

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
    auto dist_box = this->cell_matrix_inverse*point1 - this->cell_matrix_inverse*point0;

    // Extract the integer and fractional components of the distance.
    int int_x = int(dist_box.x());
    int int_y = int(dist_box.y());
    int int_z = int(dist_box.z());
    double frac_x = dist_box.x() - int_x;
    double frac_y = dist_box.y() - int_y;
    double frac_z = dist_box.z() - int_z;

    // Shift to box.

    // x
    if      (frac_x >=  0.5) frac_x -= 1.0;
    else if (frac_x <= -0.5) frac_x += 1.0;

    // y
    if      (frac_y >=  0.5) frac_y -= 1.0;
    else if (frac_y <= -0.5) frac_y += 1.0;

    // z
    if      (frac_z >=  0.5) frac_z -= 1.0;
    else if (frac_z <= -0.5) frac_z += 1.0;


    // Construct a vector from the fractional components.
    Vector frac_dist(frac_x, frac_y, frac_z);

    // Return the fractional distance vector mapped back to the triclinic
    // cell space.
    return this->cell_matrix*frac_dist;
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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        //add the delta to the coordinates of atom0
        Vector point0 = array0[i] + wrapdelta;
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            mat[j] = (array1[j] - point0);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist = qMin(mat[j].length(), mindist);
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

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group.aaBox().center(), point);
    Vector wrapped_point = point + wrapdelta;

    mat.setOuterIndex(0);

    for (int j=0; j<n; ++j)
    {
        mat[j] = (array[j] - wrapped_point);

        //store the minimum distance, the value expected to be the minimum
        //value is most efficiently placed as the second argument
        mindist = qMin(mat[j].length(), mindist);
    }

    //return the minimum distance
    return mindist;
}

/** Calculate the angle between the passed three points. This should return
    the acute angle between the points, which should lie between 0 and 180 degrees */
Angle TriclinicBox::calcAngle(const Vector &point0, const Vector &point1,
                              const Vector &point2) const
{
    Vector p0 = this->getMinimumImage(point0, point1);
    Vector p2 = this->getMinimumImage(point2, point1);

    return Vector::angle(p0, point1, p2);
}

/** Calculate the torsion angle between the passed four points. This should
    return the torsion angle measured clockwise when looking down the
    torsion from point0-point1-point2-point3. This will lie between 0 and 360
    degrees */
Angle TriclinicBox::calcDihedral(const Vector &point0, const Vector &point1,
                                 const Vector &point2, const Vector &point3) const
{
    Vector p0 = this->getMinimumImage(point0, point1);
    Vector p2 = this->getMinimumImage(point2, point1);
    Vector p3 = this->getMinimumImage(point3, point1);

    return Vector::dihedral(p0, point1, p2, p3);
}

/** Return whether or not two groups enclosed by the AABoxes 'aabox0' and
    'aabox1' are definitely beyond the cutoff distance 'dist' */
bool TriclinicBox::beyond(double dist, const AABox &aabox0, const AABox &aabox1) const
{
    Vector wrapdelta = this->wrapDelta(aabox0.center(), aabox1.center());

    return Vector::distance2(aabox0.center()+wrapdelta, aabox1.center()) >
                      SireMaths::pow_2(dist + aabox0.radius() + aabox1.radius());
}

/** Return whether or not these two groups are definitely beyond the cutoff distance. */
bool TriclinicBox::beyond(double dist, const CoordGroup &group0,
                          const CoordGroup &group1) const
{
    return this->beyond(dist, group0.aaBox(), group1.aaBox());
}

/** Return the minimum distance between the two boxes */
double TriclinicBox::minimumDistance(const AABox &box0, const AABox &box1) const
{
    // Get the distance between the minimum image of the box centers.
    Vector delta = this->wrapDelta(box0.center(), box1.center());
    delta = (box0.center() + delta) - box1.center();

    // Take absolute value of each component.
    delta = Vector(std::abs(delta.x()), std::abs(delta.y()), std::abs(delta.z()));

    // Substract box extents.
    delta -= box0.halfExtents();
    delta -= box1.halfExtents();

    // Are the boxes inside each other?
    delta = delta.max(Vector(0));

    return delta.length();
}

/** Return the minimum distance between the points in 'group0' and 'group1'.
    If this is a periodic space then this uses the minimum image convention
    (i.e. the minimum distance between the closest periodic replicas are
    used) */
double TriclinicBox::minimumDistance(const CoordGroup &group0,
                                     const CoordGroup &group1) const
{
    double mindist2(std::numeric_limits<double>::max());

    int n0 = group0.count();
    int n1 = group1.count();

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    //see if we need to wrap the coordinates...
    Vector wrapdelta = this->wrapDelta(group0.aaBox().center(), group1.aaBox().center());

    for (int i=0; i<n0; ++i)
    {
        Vector point0 = array0[i] + wrapdelta;

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two atoms
            double tmpdist = Vector::distance2(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist2 = qMin(tmpdist, mindist2);
        }
    }

    //return the minimum distance
    return sqrt(mindist2);
}

/** Return the copy of the point 'point' which is the closest minimum image
    to 'center' */
Vector TriclinicBox::getMinimumImage(const Vector &point, const Vector &center) const
{
    // Calculate the position of point in "box" space.
    auto point_box = this->cell_matrix_inverse*point;

    // Get the box shift.
    auto wrapdelta = this->wrapDelta(point, center);

    // Shift the point back to the image closest to center.
    return point + wrapdelta;
}

/** Return the closest periodic copy of 'group' to the point 'point',
    according to the minimum image convention. The effect of this is
    to move 'group' into the box which is now centered on 'point' */
CoordGroup TriclinicBox::getMinimumImage(const CoordGroup &group,
                                         const Vector &point) const
{
    Vector wrapdelta = this->wrapDelta(group.aaBox().center(), point);

    if (wrapdelta.isZero())
    {
        // Already got the minimum image.
        return group;
    }
    else
    {
        CoordGroupEditor editor = group.edit();
        editor.translate(wrapdelta);

        return editor.commit();
    }
}

/** Private function used to get the minimum image of all of the
    groups in 'groups' */
CoordGroupArray TriclinicBox::_pvt_getMinimumImage(const CoordGroupArray &groups,
                                                   const Vector &point) const
{
    int ncg = groups.count();

    const CoordGroup *group_array = groups.constData();

    if (ncg == 1)
        return CoordGroupArray( this->getMinimumImage(group_array[0], point));

    //create a new array of the right size
    QVector<CoordGroup> moved_groups(ncg);
    CoordGroup *moved_array = moved_groups.data();

    for (int i=0; i<ncg; ++i)
    {
        moved_array[i] = this->getMinimumImage(group_array[i], point);
    }

    return CoordGroupArray(moved_groups);
}

/** Return the closest periodic copy of each group in 'groups' to the
    point 'point', according to the minimum image convention.
    The effect of this is to move each 'group' into the box which is
    now centered on 'point'. If 'translate_as_one' is true,
    then this treats all groups as being part of one larger
    group, and so it translates it together. This is useful
    to get the minimum image of a molecule as a whole, rather
    than breaking the molecule across a box boundary */
CoordGroupArray TriclinicBox::getMinimumImage(const CoordGroupArray &groups,
                                              const Vector &point,
                                              bool translate_as_one) const
{
    if (translate_as_one or groups.nCoordGroups() == 1)
    {
        Vector wrapdelta = this->wrapDelta(groups.aaBox().center(), point);

        if (wrapdelta.isZero())
        {
            return groups;
        }
        else
        {
            CoordGroupArray wrapped_groups(groups);
            wrapped_groups.translate(wrapdelta);

            return wrapped_groups;
        }
    }
    else
    {
        //run through all of the groups and see if any of them need moving...
        int ncg = groups.count();

        const CoordGroup *group_array = groups.constData();

        for (int i=0; i<ncg; ++i)
        {
            const CoordGroup &group = group_array[i];

            Vector wrapdelta = this->wrapDelta(point, group.aaBox().center());

            if (not wrapdelta.isZero())
            {
                //there is at least one CoordGroup that needs moving
                // - look to translate them all!
                return this->_pvt_getMinimumImage(groups, point);
            }
        }

        //all of the CoordGroups are in the box - just return the original array
        return groups;
    }
}

/** Return the copy of the triclinic box which is the closest minimum image
    to 'center' */
AABox TriclinicBox::getMinimumImage(const AABox &aabox, const Vector &center) const
{
    Vector wrapdelta = this->wrapDelta(aabox.center(), center);

    if (wrapdelta.isZero())
        return aabox;
    else
    {
        // Get the position of the AABox center in the minimum image.
        auto min_center = this->getMinimumImage(aabox.center(), center);

        // Translate the box to this position.
        AABox ret(aabox);
        ret.translate(min_center - aabox.center());

        return ret;
    }
}

/** Return all periodic images of 'point' with respect to 'center' within
    'dist' distance of 'center' */
QVector<Vector> TriclinicBox::getImagesWithin(const Vector &point, const Vector &center,
                                              double dist) const
{
    QVector<Vector> points;

    // First, get the minimum image...
    Vector p = getMinimumImage(point, center);

    if (this->calcDist(point, center) < dist)
    {
        // The minimum image is within the distance, so lets now look at all periodic replicas...

        // We only need to look at periodic boxes that are within 'dist'...
        // This rounds to the nearest number of box lengths, e.g.
        // if dist is >= halflength.x() and < 1.5 length.x()
        // then there is only the need to go out to the first layer in the
        // x-dimension.
        int nlayers_x = int( (dist*invlength.x()) + 0.5 );
        int nlayers_y = int( (dist*invlength.y()) + 0.5 );
        int nlayers_z = int( (dist*invlength.z()) + 0.5 );

        //loop over all peridic boxes in range
        for (int i = -nlayers_x; i <= nlayers_x; ++i)
        {
            for (int j = -nlayers_y; j <= nlayers_y; ++j)
            {
                for (int k = -nlayers_z; k <= nlayers_z; ++k)
                {
                    // Get the delta value needed to translate the minimum
                    // image into the i,j,k box
                    auto delta = this->cell_matrix*Vector(i, j, k);

                    // Translate just the center of the minimum image...
                    Vector p_image = p + delta;

                    if (Vector::distance(center, p_image) < dist)
                    {
                        points.append(p_image);
                    }
                }
            }
        }
    }

    return points;
}

/** Return a list of copies of CoordGroup 'group' that are within
    'distance' of the CoordGroup 'center', translating 'group' so that
    it has the right coordinates to be around 'center'. Note that multiple
    copies of 'group' may be returned in this is a periodic space and
    there are multiple periodic replicas of 'group' within 'dist' of
    'center'. The copies of 'group' are returned together with the
    minimum distance between that periodic replica and 'center'.

    If there are no periodic replicas of 'group' that are within
    'dist' of 'center', then an empty list is returned. */
QList< tuple<double,CoordGroup> >
TriclinicBox::getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                              double dist) const
{
    if (dist > this->max_length)
        throw SireError::invalid_arg( QObject::tr(
            "You cannot use a distance (%1) that is greater than the "
            "maximum box length (%2).")
                .arg(dist).arg(this->max_length), CODELOC );

    // are there any copies within range?
    if (this->beyond(dist, group,center))
        //yep - there are no copies that are sufficiently close
        return QList< tuple<double,CoordGroup> >();

    //ok, first move 'group' into the box that has its center at the
    //same point as the center of the center group - this will give us
    //the group that is closest to us (the minimum image)
    CoordGroup minimum_image = this->getMinimumImage(group, center.aaBox().center());

    //now loop over periodic boxes, moving ever outward, trying to find
    //all copies that are within the distance

    //we can work out the maximum number of layers to go out to based on
    //the radii of the two groups, the maximum distance, and the dimensions
    //of the box
    const AABox &centerbox = center.aaBox();
    const AABox &imagebox = minimum_image.aaBox();

    double sum_of_radii_and_distance = centerbox.radius() +
                                       imagebox.radius() + dist;

    double sum_of_radii_and_distance2 = SireMaths::pow_2(sum_of_radii_and_distance);

    // this rounds to the nearest number of box lengths, e.g.
    // if sum_of_radii_and_distance is >= halflength.x() and < 1.5 length.x()
    // then there is only the need to go out to the first layer in the
    // x-dimension.
    int nlayers_x = int( (sum_of_radii_and_distance*invlength.x()) + 0.5 );
    int nlayers_y = int( (sum_of_radii_and_distance*invlength.y()) + 0.5 );
    int nlayers_z = int( (sum_of_radii_and_distance*invlength.z()) + 0.5 );

    QList< tuple<double,CoordGroup> > neargroups;

    //loop over all cubes
    for (int i = -nlayers_x; i <= nlayers_x; ++i)
    {
        for (int j = -nlayers_y; j <= nlayers_y; ++j)
        {
            for (int k = -nlayers_z; k <= nlayers_z; ++k)
            {
                // Get the delta value needed to translate the minimum
                // image into the i,j,k box
                auto delta = this->cell_matrix*Vector(i, j, k);

                //translate just the center of the minimum image...
                Vector center_of_replica = imagebox.center() + delta;

                //is the box in range?
                if ( Vector::distance2(center_of_replica,centerbox.center())
                                    <= sum_of_radii_and_distance2 )
                {
                    //yes it is! Translate the entire CoordGroup
                    CoordGroupEditor editor = minimum_image.edit();
                    editor.translate(delta);
                    CoordGroup periodic_replica = editor.commit();

                    //calculate the minimum distance... (using the cartesian space)
                    double mindist = Cartesian::minimumDistance(periodic_replica, center);

                    if (mindist <= dist)
                    {
                        neargroups.append(
                                  tuple<double,CoordGroup>(mindist,periodic_replica) );
                    }
                }
            }
        }
    }

    return neargroups;
}

/** Return a random point within the box (placing the center of the box
    is at the center 'center') */
Vector TriclinicBox::getRandomPoint(const Vector &center,
                                    const RanGenerator &generator) const
{
    return this->cell_matrix*Vector(generator.rand(-0.5, 0.5),
                                    generator.rand(-0.5, 0.5),
                                    generator.rand(-0.5, 0.5))
           + center;
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at the origin */
Vector TriclinicBox::getBoxCenter(const Vector &p) const
{
	return this->wrapDelta(Vector(0,0,0), p);
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at 'center' */
Vector TriclinicBox::getBoxCenter(const Vector &p, const Vector &center) const
{
	return this->wrapDelta(center, p);
}

const char* TriclinicBox::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TriclinicBox>() );
}
