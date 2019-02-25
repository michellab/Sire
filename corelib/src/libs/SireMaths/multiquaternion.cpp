/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2016  Christopher Woods
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

#include "multiquaternion.h"

#include <QDebug>

using namespace SireMaths;

/** Construct a null (identity) quaternion */
MultiQuaternion::MultiQuaternion()
{
    for (int i=0; i<3; i++)
        sc[i] = MultiDouble(0);

    sc[3] = MultiDouble(1);
}

/** Construct a quaternion from x,y,z,w - this normalises the values, so could be slow */
MultiQuaternion::MultiQuaternion(const MultiDouble &x, const MultiDouble &y,
                                 const MultiDouble &z, const MultiDouble &w)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
    sc[3] = w;

    MultiDouble lgth(sc[0]);
    lgth *= sc[0];
    lgth.multiplyAdd( sc[1], sc[1] );
    lgth.multiplyAdd( sc[2], sc[2] );
    lgth.multiplyAdd( sc[3], sc[3] );

    lgth = lgth.rsqrt();

    sc[0] *= lgth;
    sc[1] *= lgth;
    sc[2] *= lgth;
    sc[3] *= lgth;
}

/** Copy constructor */
MultiQuaternion::MultiQuaternion(const MultiQuaternion& p)
{
    for (int i=0; i<4; i++)
        sc[i] = p.sc[i];
}

/** Construct a quaternion which represents a rotation of 'angle' around 'axis' */
MultiQuaternion::MultiQuaternion(const MultiFloat &angle, const MultiVector &axis)
{
    //the unit quaternion can be represented by;
    // Q = cos(theta) + u*sin(theta)
    // which represents a rotation of 2*theta around the
    // vector u

    MultiFloat sintheta;
    {
        MultiFloat half_angle(angle);
        half_angle *= 0.5f;
        MultiFloat costheta;
        sincos(half_angle, sintheta, costheta);
        sc[3] = costheta;
    }

    //the vector must be normalised
    MultiVector norm = axis.normalise();

    sc[0] = sintheta * axis.x();
    sc[1] = sintheta * axis.y();
    sc[2] = sintheta * axis.z();
}

/** Construct a quaternion which represents a rotation of 'angle' around 'axis' */
MultiQuaternion::MultiQuaternion(const MultiDouble &angle, const MultiVector &axis)
{
    //the unit quaternion can be represented by;
    // Q = cos(theta) + u*sin(theta)
    // which represents a rotation of 2*theta around the
    // vector u

    MultiFloat sintheta;
    {
        MultiFloat half_angle(angle);
        half_angle *= 0.5;
        MultiFloat costheta;
        sincos(half_angle, sintheta, costheta);
        sc[3] = costheta;
    }

    //the vector must be normalised
    MultiVector norm = axis.normalise();

    sc[0] = sintheta * axis.x();
    sc[1] = sintheta * axis.y();
    sc[2] = sintheta * axis.z();
}

MultiQuaternion::~MultiQuaternion()
{}

MultiQuaternion SireMaths::operator*(const MultiQuaternion &q,
                                                      const MultiVector &p)
{
    //quaternion multiplication - p is [0, p]
    MultiDouble nw = -p.sc[0]*q.sc[0] - p.sc[1]*q.sc[1] - p.sc[2]*q.sc[2];

    //do the cross product
    MultiDouble cx = (q.sc[1]*p.z())-(q.sc[2]*p.y());
    MultiDouble cy = (q.sc[2]*p.x())-(q.sc[0]*p.z());
    MultiDouble cz = (q.sc[0]*p.y())-(q.sc[1]*p.x());

    MultiDouble nx = q.sc[3]*p.x() + cx;
    MultiDouble ny = q.sc[3]*p.y() + cy;
    MultiDouble nz = q.sc[3]*p.z() + cz;

    return MultiQuaternion(nx,ny,nz,nw);
}

MultiQuaternion SireMaths::operator*(const MultiVector &p,
                                                      const MultiQuaternion &q)
{
    return q*p;
}

/** Use this quaternion to rotate 'p' */
MultiVector MultiQuaternion::rotate(const MultiVector &p) const
{
    const MultiDouble sx2 = sc[0]*sc[0];
    const MultiDouble sy2 = sc[1]*sc[1];
    const MultiDouble sz2 = sc[2]*sc[2];

    const MultiDouble sxy = sc[0]*sc[1];
    const MultiDouble sxz = sc[0]*sc[2];
    const MultiDouble syz = sc[1]*sc[2];

    const MultiDouble swx = sc[0]*sc[3];
    const MultiDouble swy = sc[1]*sc[3];
    const MultiDouble swz = sc[2]*sc[3];

    const MultiDouble two(2.0);
    const MultiDouble half(0.5);

    return MultiVector( two*( ( half - sy2 - sz2 )*p.x()
                           + ( sxy - swz )       *p.y()
                           + ( sxz + swy )       *p.z()),

                        two*( ( sxy + swz )      *p.x()
                           + ( half - sx2 - sz2 ) *p.y()
                           + ( syz - swx )       *p.z()),

                        two*( ( sxz - swy )      *p.x()
                           + ( syz + swx )       *p.y()
                           + ( half - sx2 - sy2 ) *p.z()) );
}

/** Use the quaternion to rotate all of the points in 'p' */
QVector<MultiVector> MultiQuaternion::rotate(const QVector<MultiVector> &points) const
{
    const MultiDouble sx2 = sc[0]*sc[0];
    const MultiDouble sy2 = sc[1]*sc[1];
    const MultiDouble sz2 = sc[2]*sc[2];

    const MultiDouble sxy = sc[0]*sc[1];
    const MultiDouble sxz = sc[0]*sc[2];
    const MultiDouble syz = sc[1]*sc[2];

    const MultiDouble swx = sc[0]*sc[3];
    const MultiDouble swy = sc[1]*sc[3];
    const MultiDouble swz = sc[2]*sc[3];

    QVector<MultiVector> ret(points);
 
    const MultiDouble two(2.0);
    const MultiDouble half(0.5);
   
    for (int i=0; i<ret.count(); ++i)
    {
        MultiVector &p = ret[i];
        
        p = MultiVector( two*( ( half - sy2 - sz2 )*p.x()
                            + ( sxy - swz )      *p.y()
                            + ( sxz + swy )      *p.z()),

                         two*( ( sxy + swz )      *p.x()
                            + ( half - sx2 - sz2 ) *p.y()
                            + ( syz - swx )       *p.z()),

                         two*( ( sxz - swy )      *p.x()
                            + ( syz + swx )       *p.y()
                            + ( half - sx2 - sy2 ) *p.z()) );
    }

    return ret;
}

/** Return a quaternion that represents the identity matrix */
MultiQuaternion MultiQuaternion::identity()
{
    return MultiQuaternion();
}

/** Return a string representation of this Quaternion */
QString MultiQuaternion::toString() const
{
    QStringList parts;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        parts.append( this->at(i).toString() );
    }

    return QObject::tr("{ %1 }").arg(parts.join(", "));
}

/** Return the spherical linear interpolation (slerp) of this quaternion
    with another for 0<=lambda<=1, with this quaternion returned
    at lambda=0 and the other returned at lambda=1 */
MultiQuaternion MultiQuaternion::slerp(const MultiQuaternion &q, const MultiDouble &lam) const
{
    //slerp(q1,q2,lambda) = q1*(q1^-1*q2)^lambda
    return *this * (this->inverse()*q).pow(lam);
}

/** Return this quaternion raised to the power 'n' */
MultiQuaternion MultiQuaternion::pow(const MultiDouble &n) const
{
    // q = w + (x,y,z)
    // q = cos(theta) + u*sin(theta)
    // q^n = cos(n*theta) + u*sin(n*theta)

    MultiFloat ntheta = n*sc[3];

    MultiFloat costheta, sintheta;
    sincos(ntheta, sintheta, costheta);

    return MultiQuaternion(sc[0]*sintheta,sc[1]*sintheta,sc[2]*sintheta,costheta);
}

const MultiDouble& MultiQuaternion::x() const
{
    return sc[0];
}

const MultiDouble& MultiQuaternion::y() const
{
    return sc[1];
}

const MultiDouble& MultiQuaternion::z() const
{
    return sc[2];
}

const MultiDouble& MultiQuaternion::w() const
{
    return sc[3];
}

bool MultiQuaternion::operator==(const MultiQuaternion &p1) const
{
    return p1.sc[0] == sc[0] and p1.sc[1] == sc[1] and p1.sc[2] == sc[2] 
                  and p1.sc[3] == sc[3];
}

bool MultiQuaternion::operator!=(const MultiQuaternion &p1) const
{
    return not operator==(p1);
}

MultiQuaternion SireMaths::operator*(const MultiQuaternion &p1,
                                                      const MultiQuaternion &p2)
{
    //quaternion multiplication - if quat = [w1, v1], then

    // [w1,v1][w2,v2] = [w1w2 - v1.v2, w1v2 + w2v1 + v1xv2]

    //calculate new 'w'
    MultiDouble nw = p1.sc[3]*p2.sc[3] -
                     (p1.sc[0]*p2.sc[0] + p1.sc[1]*p2.sc[1] + p1.sc[2]*p2.sc[2]);

    //now do cross product, 'v1xv2'
    MultiDouble cx = (p1.sc[1]*p2.sc[2])-(p1.sc[2]*p2.sc[1]);
    MultiDouble cy = (p1.sc[2]*p2.sc[0])-(p1.sc[0]*p2.sc[2]);
    MultiDouble cz = (p1.sc[0]*p2.sc[1])-(p1.sc[1]*p2.sc[0]);

    MultiDouble nx = p1.sc[3]*p2.sc[0] + p2.sc[3]*p1.sc[0] + cx;
    MultiDouble ny = p1.sc[3]*p2.sc[1] + p2.sc[3]*p1.sc[1] + cy;
    MultiDouble nz = p1.sc[3]*p2.sc[2] + p2.sc[3]*p1.sc[2] + cz;

    return MultiQuaternion(nx,ny,nz,nw);
}

MultiQuaternion& MultiQuaternion::operator=(const MultiQuaternion &p)
{
    for (int i=0; i<4; i++)
        sc[i] = p.sc[i];

    return *this;
}

MultiQuaternion SireMaths::operator+(const MultiQuaternion &p1,
                                                      const MultiQuaternion &p2)
{
    return MultiQuaternion(p1.sc[0]+p2.sc[0],p1.sc[1]+p2.sc[1],
                           p1.sc[2]+p2.sc[2],p1.sc[3]+p2.sc[3]);
}

MultiQuaternion SireMaths::operator-(const MultiQuaternion &p1,
                                                      const MultiQuaternion &p2)
{
    return MultiQuaternion(p1.sc[0]-p2.sc[0],p1.sc[1]-p2.sc[1],
                           p1.sc[2]-p2.sc[2],p1.sc[3]-p2.sc[3]);
}

/** Return the inverse of the quaternion
    - since the length=1 this is the same as the conjugate */
MultiQuaternion MultiQuaternion::inverse() const
{
    MultiQuaternion ret;

    ret.sc[0] = -sc[0];
    ret.sc[1] = -sc[1];
    ret.sc[2] = -sc[2];
    ret.sc[3] = sc[3];

    return ret;
}

/** Return the conjugate of the quaternion */
MultiQuaternion MultiQuaternion::conjugate() const
{
    MultiQuaternion ret;

    ret.sc[0] = -sc[0];
    ret.sc[1] = -sc[1];
    ret.sc[2] = -sc[2];
    ret.sc[3] = sc[3];

    return ret;
}

/** Return the dot product of this with another quaternion */
MultiDouble MultiQuaternion::dot(const MultiQuaternion &q) const
{
    MultiDouble d(sc[0]);
    d *= q.sc[0];
    
    d.multiplyAdd( sc[1], q.sc[1] );
    d.multiplyAdd( sc[2], q.sc[2] );
    d.multiplyAdd( sc[3], q.sc[3] );

    return d;
}

/** Renormalise the quaternion */
void MultiQuaternion::renormalise()
{
    MultiDouble lgth(sc[0]);
    lgth *= sc[0];
    lgth.multiplyAdd( sc[1], sc[1] );
    lgth.multiplyAdd( sc[2], sc[2] );
    lgth.multiplyAdd( sc[3], sc[3] );

    lgth = lgth.rsqrt();

    sc[0] *= lgth;
    sc[1] *= lgth;
    sc[2] *= lgth;
    sc[3] *= lgth;
}

const char* MultiQuaternion::typeName()
{
    return "SireMaths::MultiQuaternion";
}

/** In-place addition operator */
MultiQuaternion& MultiQuaternion::operator+=(const MultiQuaternion &p)
{
    return this->operator=( *this + p );
}

/** In-place subtraction operator */
MultiQuaternion& MultiQuaternion::operator-=(const MultiQuaternion &p)
{
    return this->operator=( *this - p );
}

/** In-place multiplication operator */
MultiQuaternion& MultiQuaternion::operator*=(const MultiQuaternion &p)
{
    return this->operator=( *this * p );
}

/** In-place multiplication operator */
MultiQuaternion& MultiQuaternion::operator*=(const MultiVector &p)
{
    return this->operator=( *this * p );
}

/** Return the number of quaternions in this vector */
int MultiQuaternion::size()
{
    return MultiDouble::count();
}

/** Return the number of quaternions in this vector */
int MultiQuaternion::count()
{
    return MultiDouble::count();
}

/** Set the ith quaternion in this vector equal to 'val' */
void MultiQuaternion::set(int i, const Quaternion &val)
{
    sc[0].set(i, val.sc[0]);
    sc[1].quickSet(i, val.sc[1]);
    sc[2].quickSet(i, val.sc[2]);
    sc[3].quickSet(i, val.sc[3]);
}

/** Set the ith quaternion in this vector equal to 'val' */
void MultiQuaternion::quickSet(int i, const Quaternion &val)
{
    sc[0].quickSet(i, val.sc[0]);
    sc[1].quickSet(i, val.sc[1]);
    sc[2].quickSet(i, val.sc[2]);
    sc[3].quickSet(i, val.sc[3]);
}

/** Return the ith quaternion in this vector */
Quaternion MultiQuaternion::operator[](int i) const
{
    Quaternion ret;
    ret.sc[0] = sc[0].get(i);
    ret.sc[1] = sc[1].get(i);
    ret.sc[2] = sc[2].get(i);
    ret.sc[3] = sc[3].get(i);
    
    return ret;
}

/** Return the ith quaternion in this vector */
Quaternion MultiQuaternion::at(int i) const
{
    return operator[](i);
}

/** Return the ith quaternion in this vector */
Quaternion MultiQuaternion::getitem(int i) const
{
    return operator[](i);
}

/** Swap the values of the value at index idx0 in 'f0' with the value at index 'idx' in 'f1' */
void MultiQuaternion::swap(MultiQuaternion &v0, int idx0, MultiQuaternion &v1, int idx1)
{
    MultiDouble::swap(v0.sc[0], idx0, v1.sc[0], idx1);
    MultiDouble::swap(v0.sc[1], idx0, v1.sc[1], idx1);
    MultiDouble::swap(v0.sc[2], idx0, v1.sc[2], idx1);
    MultiDouble::swap(v0.sc[3], idx0, v1.sc[3], idx1);
}
