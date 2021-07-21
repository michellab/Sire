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

#include "vector.h"
#include "quaternion.h"
#include "matrix.h"
#include "nvector.h"
#include "sincos.h"

#include "SireBase/quickcopy.hpp"
#include "SireID/index.h"

#include "SireUnits/units.h"
#include "SireStream/datastream.h"

#include "SireError/errors.h"

#include <QString>
#include <QRegExp>

#include <cmath>

#include <QDebug>

using namespace SireMaths;
using namespace SireID;
using namespace SireBase;
using namespace SireUnits;
using namespace SireStream;

static const RegisterMetaType<Vector> r_vector(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const SireMaths::Vector &vec)
{
    writeHeader(ds, r_vector, 1) << vec.x() << vec.y() << vec.z();

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, SireMaths::Vector &vec)
{
    VersionID v = readHeader(ds, r_vector);

    if (v == 1)
    {
        ds >> vec.sc[0] >> vec.sc[1] >> vec.sc[2];
    }
    else
        throw version_error(v, "1", r_vector, CODELOC);

    return ds;
}

/** Construct from a tuple of three values */
Vector::Vector( const tuple<double,double,double> &pos )
{
    sc[0] = pos.get<0>();
    sc[1] = pos.get<1>();
    sc[2] = pos.get<2>();
}

/** Construct from an NVector */
Vector::Vector(const NVector &v)
{
    if (v.count() != 3)
        throw SireError::incompatible_error( QObject::tr(
                "Cannot create a 3D vector from a NVector of size %1.")
                    .arg(v.count()), CODELOC );
    
    const double *d = v.constData();
                                                    
    sc[0] = d[0];
    sc[1] = d[1];
    sc[2] = d[2];
}

/** Copy constructor */
Vector::Vector(const Vector& other)
{
    quickCopy<double>(sc, other.sc, 4);
}

/** Construct from a QString 

    \throw SireError::invalid_arg
*/
Vector::Vector(const QString &str)
{
    *this = Vector::fromString(str);
}

/** Return the inverse length squared */
double Vector::invLength2() const
{
    return double(1) / ( pow_2(sc[0]) + pow_2(sc[1]) + pow_2(sc[2]) );
}

/** Return the distance squared between two vectors */
double Vector::distance2(const Vector &v1, const Vector &v2)
{
    return pow_2(v2.sc[0]-v1.sc[0]) + pow_2(v2.sc[1]-v1.sc[1]) +
           pow_2(v2.sc[2]-v1.sc[2]);
}

/** Return the distance between two vectors */
double Vector::distance(const Vector &v1, const Vector &v2)
{
    return std::sqrt( pow_2(v2.sc[0]-v1.sc[0]) + pow_2(v2.sc[1]-v1.sc[1]) +
                      pow_2(v2.sc[2]-v1.sc[2]) );
}

/** Return the 1 / distance between two vectors */
double Vector::invDistance(const Vector &v1, const Vector &v2)
{
    return double(1.0) / std::sqrt( pow_2(v1.sc[0]-v2.sc[0]) +
                                    pow_2(v1.sc[1]-v2.sc[1]) +
                                    pow_2(v1.sc[2]-v2.sc[2]) );
}

/** Return 1 / distance2 between two vectors */
double Vector::invDistance2(const Vector &v1, const Vector &v2)
{
    return double(1.0) / ( pow_2(v1.sc[0]-v2.sc[0]) +
                           pow_2(v1.sc[1]-v2.sc[1]) +
                           pow_2(v1.sc[2]-v2.sc[2]) );
}

/** Access the elements of the array via an index operator */
double Vector::operator[](int i) const
{
    i = Index(i).map(3);
    return sc[i];
}

/** Return the size of the Vector (always 3 - unless you disagree
    with me that we should be living in a 3-dimensional space!) */
int Vector::count() const
{
    return 3;
}

/** Access elements by index */
double Vector::at(int i) const
{
    return this->operator[](i);
}

/** Access elements by index (used by python) */
double Vector::getitem(int i) const
{
    return this->operator[](i);
}

/** Return whether or not this is a zero length vector */
bool Vector::isZero() const
{
    return SireMaths::isZero(sc[0]) and SireMaths::isZero(sc[1]) 
                and SireMaths::isZero(sc[2]);
}

/** Return a normalised form of the vector */
Vector Vector::normalise() const
{
    double l = length2();

    if (SireMaths::isZero(l))
        throw SireMaths::math_error(QObject::tr(
            "Cannot normalise a zero vector, %1").arg(this->toString()),CODELOC);

    l = double(1) / sqrt(l);
    return Vector(sc[0]*l,sc[1]*l,sc[2]*l);
}

/** Return the dot product of v0 and v1 */
double Vector::dot(const Vector &v0, const Vector &v1)
{
    return (v0.sc[0]*v1.sc[0] + v0.sc[1]*v1.sc[1] + v0.sc[2]*v1.sc[2]);
}

/** Set this Vector so that it has the maximum x/y/z components out of
    this and 'other' (e.g. this->x = max(this->x(),other.x() etc.) */
void Vector::setMax(const Vector &other)
{
    for (int i=0; i<3; i++)
        sc[i] = SIRE_MAX( other.sc[i], sc[i] );
}

/** Set this Vector so that it has the minimum x/y/z components */
void Vector::setMin(const Vector &other)
{
    for (int i=0; i<3; i++)
        sc[i] = SIRE_MIN( other.sc[i], sc[i] );
}

/** Return a vector that has the maximum x/y/z components out of this
    and 'other' */
Vector Vector::max(const Vector &other) const
{
    Vector v(*this);
    v.setMax(other);
    return v;
}

/** Return a vector that has the minimum components */
Vector Vector::min(const Vector &other) const
{
    Vector v(*this);
    v.setMin(other);
    return v;
}

/** Construct a Vector from the QString representation returned by 'toString()' 

    \throw SireError::invalid_arg
*/
Vector Vector::fromString(const QString &str)
{
    QRegExp vecregexp("([0-9.-]+)\\s*,\\s*([0-9.-]+)\\s*,\\s*([0-9.-]+)");

    //use a regexp to extract the contents...
    if (vecregexp.indexIn(str) == -1)
    {
        throw SireError::invalid_arg( QObject::tr(
                "Cannot find anything that looks like a vector in the string "
                "\"%1\"").arg(str), CODELOC );
    }
        
    return Vector(vecregexp.cap(1).toFloat(),
                  vecregexp.cap(2).toFloat(),
                  vecregexp.cap(3).toFloat());
}

/** Return a QString representation of the vector */
QString Vector::toString() const
{
    return QObject::tr("( %1, %2, %3 )").arg(sc[0]).arg(sc[1]).arg(sc[2]);
}

/** Return the components via rgb (limited between 0 and 1) */
double Vector::r() const
{
    return std::max(0.0, std::min(1.0,sc[0]));
}

/** Return the components via rgb (limited between 0 and 1) */
double Vector::g() const
{
    return std::max(0.0, std::min(1.0,sc[1]));
}

/** Return the components via rgb (limited between 0 and 1) */
double Vector::b() const
{
    return std::max(0.0, std::min(1.0,sc[2]));
}

/** Set individual values of the vector */
void Vector::set(double x, double y, double z)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
}

/** Set individual values of the vector */
void Vector::set(int i, const double &v)
{
    sc[i] = v;
}

/** Set individual values of the vector */
void Vector::setX(double x)
{
    sc[0] = x;
}

/** Set individual values of the vector */
void Vector::setY(double y)
{
    sc[1] = y;
}

/** Set individual values of the vector */
void Vector::setZ(double z)
{
    sc[2] = z;
}

/** Set individual values of the vector */
void Vector::setR(double x)
{
    sc[0] = x;
}

/** Set individual values of the vector */
void Vector::setG(double y)
{
    sc[1] = y;
}

/** Set individual values of the vector */
void Vector::setB(double z)
{
    sc[2] = z;
}

/** Return the unit vector pointing in the direction of this vector */
Vector Vector::direction() const
{
    return this->normalise();
}

/** Return the length of this vector */
double Vector::magnitude() const
{
    return this->length();
}

/** Return the bearing of this vector against (0,1,0) (north) on the xy plane */
Angle Vector::bearing() const
{
    Vector t(x(),y(),0.0);

    if (t.x() == 0.0 and t.y() == 0.0)
        return Angle(0);

    if (t.y() == 0.0)
    {
       if (t.x() > 0.0)
         return 90 * degrees;
       else
         return 270 * degrees;
    }
    else if (t.y() > 0.0)  // range from -90 to 180 to 90
    {
        return (360*degrees) - Angle(std::atan(t.x()/t.y()));
    }
    else  // range from -90 to 0 to 90
    {
        return Angle(std::atan(-t.x()/t.y()));
    }
}

/** Return the bearing of this vector against 'v' on the xy plane */
Angle Vector::bearingXY(const Vector &v) const
{
    Vector px( x(), y(), 0.0);
    Vector pv(v.x(), v.y(), 0.0);
    return angle(px,pv);
}

/** Return the bearing of this vector against 'v' on the xz plane */
Angle Vector::bearingXZ(const Vector &v) const
{
    Vector px( x(), 0.0, z());
    Vector pv(v.x(), 0.0, v.z());
    return angle(px,pv);
}

/** Return the bearing of this vector against 'v' on the yz plane */
Angle Vector::bearingYZ(const Vector &v) const
{
    Vector px( 0.0, y(), z());
    Vector pv(0.0, v.y(), v.z());
    return angle(px,pv);
}

/** Return the angle between vectors 'v0' and 'v1' - this is the smallest
    angle, and will always lie between 0 and 180 degrees */
Angle Vector::angle(const Vector &v0, const Vector &v1)
{
    double d = dot(v0,v1);

    double lt = v0.length() * v1.length();

    if (SireMaths::isZero(lt))
        return Angle(0);
    else
        return Angle( std::acos(d/lt) );
}

/** Return the angle between v0-v1-v2 (treating the vectors as points in space) */
Angle Vector::angle(const Vector &v0, const Vector &v1, const Vector &v2)
{
    return angle( v0-v1, v2-v1 );
}

/** Return the dihedral angle between v0-v1-v2-v3 (treating the vectors as points) */
Angle Vector::dihedral(const Vector &v0, const Vector &v1, const Vector &v2, const Vector &v3)
{
    //     v0        v3
    //      \       /
    //      v1----v2
    //
    //    Dihedral angle is the plane of intersection between the planes
    //    formed by v0,v1,v2 and v1,v2,v3.
    //
    //    This is equivelent to the angle of intersection of the normals to those
    //    planes - thus we must calculate those normals!
    //
    //    norm0 = cross(vec(v1->v0), vec(v2->v1))
    //    norm1 = cross(vec(v2->v1), vec(v3->v2))

    //get the first plane normal...
    Vector v0_1 = v0 - v1;
    Vector v2_1 = v2 - v1;
    Vector norm0 = Vector::cross(v0_1,v2_1);

    //now get the second plane normal...
    Vector v1_2 = v1 - v2;
    Vector v3_2 = v3 - v2;
    Vector norm1 = Vector::cross(v1_2,v3_2);

    //now get the angle between the normals - this is based directly on the dot product
    //(as the normals have unit length)
    double cos_ang = Vector::dot(norm0,norm1);

    //ensure that cos_ang lies between -1 and 1
    cos_ang = SIRE_MIN(cos_ang,1.0);
    cos_ang = SIRE_MAX(cos_ang,-1.0);

    //now get the angle
    Angle ang(std::acos(cos_ang));

    //this angle will only lie between 0 and 180. To work out what
    //side of 180 this angle is we construct the plane formed by v2_1 and norm0.
    //We then calculate the angle between this new plane and that defined by norm1.
    //If this angle is < 90 degrees, then the dihedral angle is > 180 degrees

    Vector norm2 = Vector::cross(norm0,v2_1);
    cos_ang = Vector::dot(norm1,norm2);
    cos_ang = SIRE_MIN(cos_ang,1.0);
    cos_ang = SIRE_MAX(cos_ang,-1.0);

    Angle ang2(std::acos(cos_ang));

    if (ang2 < SireMaths::pi_over_two)
        return (360*degrees) - ang;
    else
        return ang;
}
    
/** Generate a vector, v0, that has distance 'dst' v0-v1, angle 'ang' v0-v1-v2,
    and dihedral 'dih' v0-v1-v2-v3 */
Vector Vector::generate(double dst, const Vector &v1, const Angle &ang, const Vector &v2,
                        const Angle &dih, const Vector &v3)
{

//     v3        v0 (generating coords of v0!)
//       \      /
//       v2---v1

//    first create a set of x/y/z orthonormal vectors, with y perpendicular
//    vec(v3-v1) and vec(v2-v1), x perpendicular to vec(v2-v1) and y, and z
//    perpendicular to x and y. Do this via cross products...

      Vector v31 = v3 - v1;
      Vector v21 = v2 - v1;

      Vector vy = Vector::cross(v31,v21);
      Vector vx = Vector::cross(v21,vy);
      Vector vz = Vector::cross(vx,vy);

//    now we have the x/y/z vectors, we can generate the new coordinates
//    from this basis set...
//    thus x/y/z in the basis set is given by
//    xbs = lgth * sin(ang) * cos(dih)
//    ybs = lgth * sin(ang) * sin(dih)
//    zbs = -lgth * cos(ang)

      double sinang, cosang;
      SireMaths::sincos(ang, &sinang, &cosang);
      double sindih, cosdih;
      SireMaths::sincos(dih, &sindih, &cosdih);

      double xbs = dst * sinang * cosdih;
      double ybs = dst * sinang * sindih;
      double zbs = -dst * cosang;

//    Then we map the coordinates in this basis set to our cartesian coordinates
//    via...
//    x = xbs*vx(1) + ybs*vy(1) + zbs*vz(1)
//    y = xbs*vx(2) + ybs*vy(2) + zbs*vz(2)
//    z = xbs*vx(3) + ybs*vy(3) + zbs*vz(3)
//
//    These coordinates are based at the origin - they need to be based from
//    the coordinates of the bond atom by adding on v1.
//
//    (we combine the last two steps together for speed)

      double nx = vx.x()*xbs + vy.x()*ybs + vz.x()*zbs + v1.x();
      double ny = vx.y()*xbs + vy.y()*ybs + vz.y()*zbs + v1.y();
      double nz = vx.z()*xbs + vy.z()*ybs + vz.z()*zbs + v1.z();

      return Vector(nx,ny,nz);
}

 /** Return the cross product of v0 and v1 
     N.B. This function returns the normalised cross product.
 */
Vector Vector::cross(const Vector &v0, const Vector &v1)
{
    double nx = v0.sc[1]*v1.sc[2] - v0.sc[2]*v1.sc[1];
    double ny = v0.sc[2]*v1.sc[0] - v0.sc[0]*v1.sc[2];
    double nz = v0.sc[0]*v1.sc[1] - v0.sc[1]*v1.sc[0];

    Vector vec(nx,ny,nz);
    double length = std::sqrt(nx*nx + ny*ny + nz*nz);
    
    if (SireMaths::isZero(length))
    {
        //these two vectors are parallel - we need to choose just
        //one perpendicular vector
        if (v0.isZero())
        {
            if (v1.isZero())
                //both are null vectors. Just return a unit vector along the x axis
                return Vector(1.0, 0.0, 0.0);
            else
            {
                //do this by creating a copy of v0 with two elements swapped
                if ( SireMaths::isZero(v1.x()) )
                    return Vector(0,v1.z(),-v1.y()).normalise();
                else
                    return Vector(v1.y(),-v1.x(),0).normalise();
            }
        }
        else
        {
            //do this by creating a copy of v0 with two elements swapped
            if ( SireMaths::isZero(v0.x()) )
                return Vector(0,v0.z(),-v0.y()).normalise();
            else
                return Vector(v0.y(),-v0.x(),0).normalise();
        }
    }
    else
    {
        double invlength = double(1)/length;
        Vector normal( vec.x()*invlength, vec.y()*invlength, vec.z()*invlength );
        double checklength = normal.x()*normal.x() + normal.y()*normal.y() + normal.z()*normal.z();
        
        if (checklength > 1.1 or checklength < 0.9)
        {
            //something went wrong with normalisation
            qDebug() << "WEIRD ERROR WITH NORMALISING VECTOR" << vec.toString() << checklength;
            qDebug() << normal.toString();
            qDebug() << "AUTOMATICALLY FIXING THE ERROR...";

            //these two vectors are parallel - we need to choose just
            //one perpendicular vector
            if (v0.isZero())
            {
                if (v1.isZero())
                    //both are null vectors. Just return a unit vector along the x axis
                    normal = Vector(1.0, 0.0, 0.0);
                else
                {
                    //do this by creating a copy of v0 with two elements swapped
                    if ( SireMaths::isZero(v1.x()) )
                        normal = Vector(v1.x(),v1.z(),v1.y()).normalise();
                    else
                        normal = Vector(v1.y(),v1.x(),v1.z()).normalise();
                }
            }
            else
            {
                //do this by creating a copy of v0 with two elements swapped
                if ( SireMaths::isZero(v0.x()) )
                    normal = Vector(v0.x(),v0.z(),v0.y()).normalise();
                else
                    normal =Vector(v0.y(),v0.x(),v0.z()).normalise();
            }

            qDebug() << "Fixed cross product of vectors" << v0.toString() << ":" << v1.toString();
            qDebug() << "equals" << normal.toString();
        }
        
        return normal;
    }
}

/** The actual cross productor of vector v0 and v1.
 */
Vector Vector::realCross(const Vector& v0, const Vector& v1)
{
    double nx = v0.y()*v1.z() - v0.z()*v1.y();
    double ny = v0.z()*v1.x() - v0.x()*v1.z();
    double nz = v0.x()*v1.y() - v0.y()*v1.x();

    return Vector(nx,ny,nz);
}

/** Return the manhattan length of the vector */
double Vector::manhattanLength() const
{
    return std::abs(sc[0])+std::abs(sc[1])+std::abs(sc[2]);
}

/** Return the metric tensor of a vector, i.e.
          
    | y*y + z*z,    -x*y    -x*z      |
    |    -y*x,   x*x + z*z  -y*z      |
    |    -z*x       -z*y    x*x + y*y |

*/
Matrix Vector::metricTensor() const
{
    double x2 = sc[0]*sc[0];
    double y2 = sc[1]*sc[1];
    double z2 = sc[2]*sc[2];

    double xy = sc[0]*sc[1];
    double xz = sc[0]*sc[2];
    double yz = sc[1]*sc[2];

    return Matrix( y2 + z2, -xy, -xz,
                   -xy, x2 + z2, -yz,
                   -xz, -yz, x2 + y2 );
}

/** Return the multiple of this vector with the matrix 'm' */
const Vector SireMaths::operator*(const Matrix &m, const Vector &p)
{
    return Vector(m.xx()*p.sc[0] + m.xy()*p.sc[1] + m.xz()*p.sc[2],
                  m.yx()*p.sc[0] + m.yy()*p.sc[1] + m.yz()*p.sc[2],
                  m.zx()*p.sc[0] + m.zy()*p.sc[1] + m.zz()*p.sc[2]);
}

/** Increment, decrement, negate etc. */
const Vector& Vector::operator/=(const double &val)
{
    if ( SireMaths::isZero(val) )
        throw SireMaths::math_error(QObject::tr(
            "Cannot divide a vector by zero! %1 / 0 is a error!").arg(this->toString()),CODELOC);

    for (int i=0; i<3; i++)
        sc[i] /= val;

    return *this;
}

/** Increment, decrement, negate etc. */
const Vector SireMaths::operator/(const Vector &p1, double c)
{
    if (isZero(c))
        throw SireMaths::math_error(QObject::tr(
            "Cannot divide a vector by zero! %1 / 0 is a error!").arg(p1.toString()),CODELOC);

    return Vector(p1.sc[0]/c, p1.sc[1]/c, p1.sc[2]/c);
}

uint get_hash(double val)
{
    return qHash( (quint64)(val) );
}

uint qHash(const SireMaths::Vector &vec)
{
    return get_hash(vec.x()) + get_hash(vec.y()) + get_hash(vec.z());
}

const char* Vector::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Vector>() );
}
