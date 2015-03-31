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

#include "quaternion.h"
#include "vector.h"
#include "matrix.h"

#include <QString>
#include <QRegExp>
#include <cmath>

#include <QDebug>

#include "SireStream/datastream.h"

using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<Quaternion> r_quat(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Quaternion &quat)
{
    writeHeader(ds, r_quat, 1)  << quat.sc[0] << quat.sc[1]
                                << quat.sc[2] << quat.sc[3];

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Quaternion &quat)
{
    VersionID v = readHeader(ds, r_quat);

    if (v == 1)
    {
        ds >> quat.sc[0] >> quat.sc[1] >> quat.sc[2] >> quat.sc[3];
    }
    else
        throw version_error(v, "1", r_quat, CODELOC);

    return ds;
}

/** Construct a null (identity) quaternion */
Quaternion::Quaternion()
{
    for (int i=0; i<3; i++)
        sc[i] = 0.0;

    sc[3] = 1.0;
}

/** Construct a quaternion from x,y,z,w - this normalises the values, so could be slow */
Quaternion::Quaternion(double x, double y, double z, double w)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
    sc[3] = w;

    double lgth2 = sc[0]*sc[0] + sc[1]*sc[1] + sc[2]*sc[2] + sc[3]*sc[3];

    if (lgth2 != double(1.0))
    {
        double invlgth = 1.0 / std::sqrt(lgth2);
        sc[0] *= invlgth;
        sc[1] *= invlgth;
        sc[2] *= invlgth;
        sc[3] *= invlgth;
    }
}

/** Copy constructor */
Quaternion::Quaternion(const Quaternion& p)
{
    for (int i=0; i<4; i++)
        sc[i] = p.sc[i];
}

/** Construct a quaternion which represents a rotation of 'angle' around 'axis' */
Quaternion::Quaternion(SireUnits::Dimension::Angle angle, const Vector &axis)
{
    //the unit quaternion can be represented by;
    // Q = cos(theta) + u*sin(theta)
    // which represents a rotation of 2*theta around the
    // vector u

    double ang = double(angle) * 0.5;
    double costheta = cos(ang);
    double sintheta = sin(ang);

    sc[3] = costheta;

    //the vector must be normalised
    double lngth2 = axis.length2();

    if (lngth2 != 1)
    {
        Vector a2 = axis.normalise();

        sc[0] = sintheta * a2.x();
        sc[1] = sintheta * a2.y();
        sc[2] = sintheta * a2.z();
    }
    else
    {
        sc[0] = sintheta * axis.x();
        sc[1] = sintheta * axis.y();
        sc[2] = sintheta * axis.z();
    }

    renormalise();
}

/** Construct from a Matrix */
Quaternion::Quaternion(const Matrix &m)
{
    this->fromMatrix(m);
}

Quaternion::~Quaternion()
{}

bool Quaternion::isIdentity() const
{
    return sc[3] == 1.0;
}

const Quaternion SIREMATHS_EXPORT SireMaths::operator*(const Quaternion &q, const Vector &p)
{
    //quaternion multiplication - p is [0, p]
    double nw = -p.sc[0]*q.sc[0] - p.sc[1]*q.sc[1] - p.sc[2]*q.sc[2];

    //do the cross product
    double cx = (q.sc[1]*p.z())-(q.sc[2]*p.y());
    double cy = (q.sc[2]*p.x())-(q.sc[0]*p.z());
    double cz = (q.sc[0]*p.y())-(q.sc[1]*p.x());

    double nx = q.sc[3]*p.x() + cx;
    double ny = q.sc[3]*p.y() + cy;
    double nz = q.sc[3]*p.z() + cz;

    return Quaternion(nx,ny,nz,nw);
}

const Quaternion SIREMATHS_EXPORT SireMaths::operator*(const Vector &p, const Quaternion &q)
{
    return q*p;
}

/** Convert into a matrix */
Matrix Quaternion::toMatrix() const
{
    /* Thanks to http://www.flipcode.com/documents/matrfaq.html#Q54
    
        |       2     2                                |
        | 1 - 2Y  - 2Z    2XY - 2ZW      2XZ + 2YW     |
        |                                              |
        |                       2     2                |
    M = | 2XY + 2ZW       1 - 2X  - 2Z   2YZ - 2XW     |
        |                                              |
        |                                      2     2 |
        | 2XZ - 2YW       2YZ + 2XW      1 - 2X  - 2Y  |
        |                                              |
    */

    return 
     Matrix( 1.0-(2.0*sc[1]*sc[1])-(2.0*sc[2]*sc[2]), (2.0*sc[0]*sc[1])-(2.0*sc[3]*sc[2]),
             (2.0*sc[0]*sc[2])+(2.0*sc[1]*sc[3]),
             (2.0*sc[0]*sc[1])+(2.0*sc[3]*sc[2]), 1.0-(2.0*sc[0]*sc[0])-(2.0*sc[2]*sc[2]),
             (2.0*sc[1]*sc[2])-(2.0*sc[3]*sc[0]),
             (2.0*sc[0]*sc[2])-(2.0*sc[1]*sc[3]), (2.0*sc[1]*sc[2])+(2.0*sc[0]*sc[3]),
             1.0-(2.0*sc[0]*sc[0])-(2.0*sc[1]*sc[1]) );
}

/** Use this quaternion to rotate 'p' */
Vector Quaternion::rotate(const Vector &p) const
{
    double sx2 = sc[0]*sc[0];
    double sy2 = sc[1]*sc[1];
    double sz2 = sc[2]*sc[2];

    double sxy = sc[0]*sc[1];
    double sxz = sc[0]*sc[2];
    double syz = sc[1]*sc[2];

    double swx = sc[0]*sc[3];
    double swy = sc[1]*sc[3];
    double swz = sc[2]*sc[3];

    return Vector( 2.0*( ( 0.5 - sy2 - sz2 )*p.x()
                       + ( sxy - swz )      *p.y()
                       + ( sxz + swy )      *p.z()),

                   2.0*( ( sxy + swz )      *p.x()
                       + ( 0.5 - sx2 - sz2 )*p.y()
                       + ( syz - swx )      *p.z()),

                   2.0*( ( sxz - swy )      *p.x()
                       + ( syz + swx )      *p.y()
                       + ( 0.5 - sx2 - sy2 )*p.z()) );
}

/** Use the quaternion to rotate all of the points in 'p' */
QVector<Vector> Quaternion::rotate(const QVector<Vector> &points) const
{
    if (this->isIdentity())
        return points;

    const double sx2 = sc[0]*sc[0];
    const double sy2 = sc[1]*sc[1];
    const double sz2 = sc[2]*sc[2];

    const double sxy = sc[0]*sc[1];
    const double sxz = sc[0]*sc[2];
    const double syz = sc[1]*sc[2];

    const double swx = sc[0]*sc[3];
    const double swy = sc[1]*sc[3];
    const double swz = sc[2]*sc[3];

    QVector<Vector> ret(points);
    
    for (int i=0; i<ret.count(); ++i)
    {
        Vector &p = ret[i];
        
        p = Vector( 2.0*( ( 0.5 - sy2 - sz2 ) *p.x()
                         + ( sxy - swz )      *p.y()
                         + ( sxz + swy )      *p.z()),

                    2.0*( ( sxy + swz )       *p.x()
                        + ( 0.5 - sx2 - sz2 ) *p.y()
                        + ( syz - swx )       *p.z()),

                    2.0*( ( sxz - swy )       *p.x()
                        + ( syz + swx )       *p.y()
                        + ( 0.5 - sx2 - sy2 ) *p.z()) );
    }

    return ret;
}

/** Return a quaternion that represents the identity matrix */
Quaternion Quaternion::identity()
{
    return Quaternion();
}

/** Get from a matrix */
void Quaternion::fromMatrix(const Matrix &m)
{
    if (not SireMaths::areEqual(m.determinant(), 1.0))
        throw SireMaths::domain_error( QObject::tr(
                "You can only convert a rotation matrix to a quaternion. Rotation matrices "
                "have a determinant of 1.0, while the matrix you passed has a determinant of "
                "%1\n%2").arg(m.determinant()).arg(m.toString()), CODELOC );

  /* Thanks to http://www.flipcode.com/documents/matrfaq.html#Q54
  
  A rotation may be converted back to a quaternion through the use of
  the following algorithm:

  The process is performed in the following stages, which are as follows:

    Calculate the trace of the matrix T from the equation:

      T = mat[0] + mat[5] + mat[10] + 1

   */
    const double trace = m.xx() + m.yy() + m.zz() + 1;

    /*
    If the trace of the matrix is greater than zero, then
    perform an "instant" calculation.

      S = 0.5 / sqrt(T)

      W = 0.25 / S                           Mapping from Sire to mat[X]
                                             xx xy xz   0 1 2
      X = ( mat[9] - mat[6] ) * S            yx yy yz   4 5 6
                                             zx zy zz   8 9 10
      Y = ( mat[2] - mat[8] ) * S

      Z = ( mat[4] - mat[1] ) * S
    */
    if (trace > 0.00001)  // use 0.00001 to avoid numerical instability near 0
    {
        double s = 0.5 / sqrt(trace);
        sc[3] = 0.25 / s;

        sc[0] = (m.zy()-m.yz())*s;
        sc[1] = (m.xz()-m.zx())*s;
        sc[2] = (m.yx()-m.xy())*s;
    }
    else
    {
        /*
        If the trace of the matrix is less than or equal to zero
        then identify which major diagonal element has the greatest
        value. */

        if ( (m.zz() >= m.yy()) and (m.zz() >= m.xx()) )
        {
            /* Column 2:
            S  = sqrt( 1.0 + mr[10] - mr[0] - mr[5] ) * 2;

            Qx = (mr[2] + mr[8] ) / S;
            Qy = (mr[6] + mr[9] ) / S;
            Qz = 0.5 / S;
            Qw = (mr[1] + mr[4] ) / S; */

            double s = sqrt( 1 + m.zz() - m.xx() - m.yy() ) * 2;

            sc[0] = (m.xz() + m.zx()) / s;
            sc[1] = (m.zy() + m.yz()) / s;
            sc[2] = 0.5 / s;
            sc[3] = (m.yx() + m.xy()) / s;
        }
        else if (m.yy() >= m.xx())
        {
           /* Column 1:
            S  = sqrt( 1.0 + mr[5] - mr[0] - mr[10] ) * 2;

            Qx = (mr[1] + mr[4] ) / S;
            Qy = 0.5 / S;
            Qz = (mr[6] + mr[9] ) / S;
            Qw = (mr[2] + mr[8] ) / S; */

            double s = sqrt( 1.0 + m.yy() - m.zz() - m.xx() ) * 2;

            sc[0] = (m.yx() + m.xy()) / s;
            sc[1] = 0.5 / s;
            sc[2] = (m.zy() + m.yz()) / s;
            sc[3] = (m.xz() + m.zx()) / s;
        }
        else
        {
            /* Column 0:
            S  = sqrt( 1.0 + mr[0] - mr[5] - mr[10] ) * 2;

            Qx = 0.5 / S;
            Qy = (mr[1] + mr[4] ) / S;
            Qz = (mr[2] + mr[8] ) / S;
            Qw = (mr[6] + mr[9] ) / S; */

            double s = sqrt( 1.0 + m.xx() - m.yy() - m.zz() ) * 2;
            sc[0] = 0.5 / s;
            sc[1] = (m.yx() + m.xy()) / s;
            sc[2] = (m.xz() + m.zx()) / s;
            sc[3] = (m.zy() + m.yz()) / s;
        }
    }

    /*
     The quaternion is then defined as:

       Q = | Qx Qy Qz Qw |
    */

    renormalise();
}

/** Return a string representation of this Quaternion */
QString Quaternion::toString() const
{
    QString s;
    return s.sprintf("(%8.3f, %8.3f, %8.3f, %8.3f)",sc[0],sc[1],sc[2],sc[3]);
}

QRegExp quatregexp("([0-9.-]+),\\s{0,}([0-9.-]+),\\s{0,}([0-9.-]+),\\s{0,}([0-9.-]+)");

/** Construct a Vector from the QString representation returned by 'toString()' */
Quaternion Quaternion::fromString(const QString &str)
{
    //use a regexp to extract the contents...
    if (quatregexp.indexIn(str) != -1)
    {
        Quaternion q(quatregexp.cap(1).toFloat(),
                    quatregexp.cap(2).toFloat(),
                    quatregexp.cap(3).toFloat(),
                    quatregexp.cap(4).toFloat());
        q.renormalise();
        return q;
    }
    else
       return Quaternion();
}

/** Return the spherical linear interpolation (slerp) of this quaternion
    with another for 0<=lambda<=1, with this quaternion returned
    at lambda=0 and the other returned at lambda=1 */
Quaternion Quaternion::slerp(const Quaternion &q, double lambda) const
{
    //slerp(q1,q2,lambda) = q1*(q1^-1*q2)^lambda

    if (lambda <= 0.0)
      return *this;
    else if (lambda >= 1.0)
      return q;
    else
    {
        return *this * (this->inverse()*q).pow(lambda);
    }
}

/** Return this quaternion raised to the power 'n' */
Quaternion Quaternion::pow(double n) const
{
    // q = w + (x,y,z)
    // q = cos(theta) + u*sin(theta)
    // q^n = cos(n*theta) + u*sin(n*theta)

    double ntheta = n*sc[3];
    double costheta = cos(ntheta);
    double sintheta = sin(ntheta);
    return Quaternion(sc[0]*sintheta,sc[1]*sintheta,sc[2]*sintheta,costheta);
}

double Quaternion::x() const
{
    return sc[0];
}

double Quaternion::y() const
{
    return sc[1];
}

double Quaternion::z() const
{
    return sc[2];
}

double Quaternion::w() const
{
    return sc[3];
}

bool Quaternion::operator==(const Quaternion &p1) const
{
    return p1.sc[0] == sc[0] and p1.sc[1] == sc[1] and p1.sc[2] == sc[2] 
                  and p1.sc[3] == sc[3];
}

bool Quaternion::operator!=(const Quaternion &p1) const
{
    return p1.sc[0] != sc[0] or p1.sc[1] != sc[1] or p1.sc[2] != sc[2] 
                  or p1.sc[3] != sc[3];
}

const Quaternion SIREMATHS_EXPORT SireMaths::operator*(const Quaternion &p1,
                                                       const Quaternion &p2)
{
    //quaternion multiplication - if quat = [w1, v1], then

    // [w1,v1][w2,v2] = [w1w2 - v1.v2, w1v2 + w2v1 + v1xv2]

    //calculate new 'w'
    double nw = p1.sc[3]*p2.sc[3] - (p1.sc[0]*p2.sc[0] + p1.sc[1]*p2.sc[1] + p1.sc[2]*p2.sc[2]);

    //now do cross product, 'v1xv2'
    double cx = (p1.sc[1]*p2.sc[2])-(p1.sc[2]*p2.sc[1]);
    double cy = (p1.sc[2]*p2.sc[0])-(p1.sc[0]*p2.sc[2]);
    double cz = (p1.sc[0]*p2.sc[1])-(p1.sc[1]*p2.sc[0]);

    double nx = p1.sc[3]*p2.sc[0] + p2.sc[3]*p1.sc[0] + cx;
    double ny = p1.sc[3]*p2.sc[1] + p2.sc[3]*p1.sc[1] + cy;
    double nz = p1.sc[3]*p2.sc[2] + p2.sc[3]*p1.sc[2] + cz;

    return Quaternion(nx,ny,nz,nw);
}

Quaternion& Quaternion::operator=(const Quaternion &p)
{
    for (int i=0; i<4; i++)
        sc[i] = p.sc[i];

    return *this;
}

const Quaternion SIREMATHS_EXPORT SireMaths::operator+(const Quaternion &p1,
                                                       const Quaternion &p2)
{
    return Quaternion(p1.sc[0]+p2.sc[0],p1.sc[1]+p2.sc[1],p1.sc[2]+p2.sc[2],p1.sc[3]+p2.sc[3]);
}

const Quaternion SIREMATHS_EXPORT SireMaths::operator-(const Quaternion &p1,
                                                       const Quaternion &p2)
{
    return Quaternion(p1.sc[0]-p2.sc[0],p1.sc[1]-p2.sc[1],p1.sc[2]-p2.sc[2],p1.sc[3]-p2.sc[3]);
}

/** Return the inverse of the quaternion
    - since the length=1 this is the same as the conjugate */
Quaternion Quaternion::inverse() const
{
    return Quaternion(-sc[0],-sc[1],-sc[2],sc[3]);
}

/** Return the conjugate of the quaternion */
Quaternion Quaternion::conjugate() const
{
    return Quaternion(-sc[0],-sc[1],-sc[2],sc[3]);
}

/** Return the dot product of this with another quaternion */
double Quaternion::dot(const Quaternion &q) const
{
    return sc[0]*q.sc[0] + sc[1]*q.sc[1] + sc[2]*q.sc[2] + sc[3]*q.sc[3];
}

/** Renormalise the quaternion */
void Quaternion::renormalise()
{
    double l = invSqrt(sc[0]*sc[0]+sc[1]*sc[1]+sc[2]*sc[2]+sc[3]*sc[3]);
    sc[0] *= l;
    sc[1] *= l;
    sc[2] *= l;
    sc[3] *= l;
}

const char* Quaternion::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Quaternion>() );
}

/** In-place addition operator */
Quaternion& Quaternion::operator+=(const Quaternion &p)
{
    return this->operator=( *this + p );
}

/** In-place subtraction operator */
Quaternion& Quaternion::operator-=(const Quaternion &p)
{
    return this->operator=( *this - p );
}

/** In-place multiplication operator */
Quaternion& Quaternion::operator*=(const Quaternion &p)
{
    return this->operator=( *this * p );
}

/** In-place multiplication operator */
Quaternion& Quaternion::operator*=(const Vector &p)
{
    return this->operator=( *this * p );
}
