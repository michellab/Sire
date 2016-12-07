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

#include "multivector.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include <QDebug>

using namespace SireMaths;
using namespace SireID;

/** Return the number of vectors in this MultiVector */
int MultiVector::size()
{
    return MultiDouble::size();
}

/** Return the number of vectors in this MultiVector */
int MultiVector::count()
{
    return MultiDouble::size();
}

/** Construct from an array of Vectors. This array cannot be greater
    than MultiDouble::size() */
MultiVector::MultiVector(const Vector *array, int sz)
{
    if (sz > MultiDouble::size())
        throw SireError::unsupported( QObject::tr(
                "Cannot fit an array of size %1 in this MultiVector, as it is only "
                "capable of holding %2 values...").arg(sz).arg(MultiDouble::size()), CODELOC );

    for (int i=0; i<sz; ++i)
    {
        sc[0].set(i, array[i].x());
        sc[1].set(i, array[i].y());
        sc[2].set(i, array[i].z());
    }
}

/** Construct from an array of Vectors. This array cannot be greater
    than MultiDouble::size() */
MultiVector::MultiVector(const QVector<Vector> &array)
{
    const int sz = array.count();

    if (sz > MultiDouble::size())
        throw SireError::unsupported( QObject::tr(
                "Cannot fit an array of size %1 in this MultiVector, as it is only "
                "capable of holding %2 values...").arg(sz).arg(MultiDouble::size()), CODELOC );

    for (int i=0; i<sz; ++i)
    {
        sc[0].set(i, array.constData()[i].x());
        sc[1].set(i, array.constData()[i].y());
        sc[2].set(i, array.constData()[i].z());
    }
}

/** Copy constructor */
MultiVector::MultiVector(const MultiVector& other)
{
    sc[0] = other.sc[0];
    sc[1] = other.sc[1];
    sc[2] = other.sc[2];
}

/** Destructor */
MultiVector::~MultiVector()
{}

/** Convert the passed array of vectors into an array of MultiVectors */
QVector<MultiVector> MultiVector::fromArray(const Vector *array, int size)
{
    if (size == 0)
        return QVector<MultiVector>();
    
    int nvecs = size / count();
    int nremain = size % count();
    
    QVector<MultiVector> marray(nvecs + ( (nremain > 0) ? 1 : 0 ));
    
    MultiVector *ma = marray.data();

    int idx = 0;

    for (int i=0; i<nvecs; ++i)
    {
        ma[i] = MultiVector(array+idx, count());
        idx += count();
    }

    if (nremain > 0)
    {
        ma[marray.count()-1] = MultiVector(array+idx, nremain);
    }
    
    return marray;
}

/** Convert the passed array of vectors into an array of MultiVectors */
QVector<MultiVector> MultiVector::fromArray(const QVector<Vector> &array)
{
    return fromArray(array.constData(), array.count());
}

/** Return the inverse length squared */
MultiDouble MultiVector::invLength2() const
{
    return length2().reciprocal();
}

/** Return the distance squared between two vectors */
MultiDouble MultiVector::distance2(const MultiVector &v1, const MultiVector &v2)
{
    MultiDouble del = v1.sc[0] - v2.sc[0];
    MultiDouble dist2( del * del );
    
    del = v1.sc[1] - v2.sc[1];
    dist2.multiplyAdd( del, del );
    
    del = v1.sc[2] - v2.sc[2];
    dist2.multiplyAdd( del, del );

    return dist2;
}

/** Return the distance between two vectors */
MultiDouble MultiVector::distance(const MultiVector &v1, const MultiVector &v2)
{
    return distance2(v1, v2).sqrt();
}

/** Return the 1 / distance between two vectors */
MultiDouble MultiVector::invDistance(const MultiVector &v1, const MultiVector &v2)
{
    return distance(v1,v2).reciprocal();
}

/** Return 1 / distance2 between two vectors */
MultiDouble MultiVector::invDistance2(const MultiVector &v1, const MultiVector &v2)
{
    return distance2(v1,v2).reciprocal();
}

/** Access the 'ith' vector in the MultiVector */
Vector MultiVector::operator[](int i) const
{
    return Vector( sc[0].get(i), sc[1].get(i), sc[2].get(i) );
}

/** Access the 'ith' vector in the MultiVector */
Vector MultiVector::at(int i) const
{
    return this->operator[](i);
}

/** Access the 'ith' vector in the MultiVector */
Vector MultiVector::getitem(int i) const
{
    return this->operator[](i);
}

/** Return a normalised form of the vector */
MultiVector MultiVector::normalise() const
{
    MultiDouble l = length();

    MultiDouble mask = l.compareEqual( MultiDouble(0) );

    l = mask.logicalAnd( l.reciprocal() );
    
    return MultiVector(sc[0]*l,sc[1]*l,sc[2]*l);
}

/** Return the dot product of v0 and v1 */
MultiDouble MultiVector::dot(const MultiVector &v0, const MultiVector &v1)
{
    MultiDouble result(v0.sc[0]);
    result *= v1.sc[0];
    result.multiplyAdd( v0.sc[1], v1.sc[1] );
    result.multiplyAdd( v0.sc[2], v1.sc[2] );

    return result;
}

/** Set this Vector so that it has the maximum x/y/z components out of
    this and 'other' (e.g. this->x = max(this->x(),other.x() etc.) */
void MultiVector::setMax(const MultiVector &other)
{
    for (int i=0; i<3; i++)
        sc[i] = sc[i].max(other.sc[i]);
}

/** Set this Vector so that it has the minimum x/y/z components */
void MultiVector::setMin(const MultiVector &other)
{
    for (int i=0; i<3; i++)
        sc[i] = sc[i].min(other.sc[i]);
}

/** Return a vector that has the maximum x/y/z components out of this
    and 'other' */
MultiVector MultiVector::max(const MultiVector &other) const
{
    MultiVector v(*this);
    v.setMax(other);
    return v;
}

/** Return a vector that has the minimum components */
MultiVector MultiVector::min(const MultiVector &other) const
{
    MultiVector v(*this);
    v.setMin(other);
    return v;
}

/** Return a QString representation of the vector */
QString MultiVector::toString() const
{
    QStringList parts;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        parts.append( this->at(i).toString() );
    }

    return QObject::tr("{ %1 }").arg(parts.join(", "));
}

/** Return the components via rgb (limited between 0 and 1) */
MultiDouble MultiVector::r() const
{
    return sc[0].max( MultiDouble(0) ).min( MultiDouble(1) );
}

/** Return the components via rgb (limited between 0 and 1) */
MultiDouble MultiVector::g() const
{
    return sc[1].max( MultiDouble(0) ).min( MultiDouble(1) );
}

/** Return the components via rgb (limited between 0 and 1) */
MultiDouble MultiVector::b() const
{
    return sc[2].max( MultiDouble(0) ).min( MultiDouble(1) );
}

/** Set individual values of the vector */
void MultiVector::set(const MultiDouble &x, const MultiDouble &y, const MultiDouble &z)
{
    sc[0] = x;
    sc[1] = y;
    sc[2] = z;
}

/** Set individual values of the vector */
void MultiVector::set(int i, const Vector &v)
{
    i = Index(i).map(MultiDouble::count());
    sc[0].quickSet(i, v.x());
    sc[1].quickSet(i, v.y());
    sc[2].quickSet(i, v.z());
}

/** Quickly set the values of the vector, without checking the index is valid */
void MultiVector::quickSet(int i, const Vector &v)
{
    sc[0].quickSet(i, v.x());
    sc[1].quickSet(i, v.y());
    sc[2].quickSet(i, v.z());
}

/** Set individual values of the vector */
void MultiVector::setX(const MultiDouble &x)
{
    sc[0] = x;
}

/** Set individual values of the vector */
void MultiVector::setY(const MultiDouble &y)
{
    sc[1] = y;
}

/** Set individual values of the vector */
void MultiVector::setZ(const MultiDouble &z)
{
    sc[2] = z;
}

/** Set individual values of the vector */
void MultiVector::setR(const MultiDouble &r)
{
    sc[0] = r.max(MultiDouble(0)).min(MultiDouble(1));
}

/** Set individual values of the vector */
void MultiVector::setG(const MultiDouble &g)
{
    sc[1] = g.max(MultiDouble(0)).min(MultiDouble(1));
}

/** Set individual values of the vector */
void MultiVector::setB(const MultiDouble &b)
{
    sc[2] = b.max(MultiDouble(0)).min(MultiDouble(1));
}

/** Return the unit vector pointing in the direction of this vector */
MultiVector MultiVector::direction() const
{
    return this->normalise();
}

/** Return the length of this vector */
MultiDouble MultiVector::magnitude() const
{
    return this->length();
}

/** Return the bearing of this vector against (0,1,0) (north) on the xy plane */
MultiDouble MultiVector::bearing() const
{
    MultiDouble result;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        result.set(i, this->at(i).bearing());
    }

    return result;
}

/** Return the bearing of this vector against 'v' on the xy plane */
MultiDouble MultiVector::bearingXY(const MultiVector &v) const
{
    MultiVector px( x(), y(), 0.0);
    MultiVector pv(v.x(), v.y(), 0.0);
    return angle(px,pv);
}

/** Return the bearing of this vector against 'v' on the xz plane */
MultiDouble MultiVector::bearingXZ(const MultiVector &v) const
{
    MultiVector px( x(), 0.0, z());
    MultiVector pv(v.x(), 0.0, v.z());
    return angle(px,pv);
}

/** Return the bearing of this vector against 'v' on the yz plane */
MultiDouble MultiVector::bearingYZ(const MultiVector &v) const
{
    MultiVector px( 0.0, y(), z());
    MultiVector pv(0.0, v.y(), v.z());
    return angle(px,pv);
}

/** Return the angle between vectors 'v0' and 'v1' - this is the smallest
    angle, and will always lie between 0 and 180 degrees */
MultiDouble MultiVector::angle(const MultiVector &v0, const MultiVector &v1)
{
    MultiDouble result;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        result.set(i, Vector::angle(v0.at(i),v1.at(i)));
    }

    return result;
}

/** Return the angle between v0-v1-v2 (treating the vectors as points in space) */
MultiDouble MultiVector::angle(const MultiVector &v0, const MultiVector &v1, const MultiVector &v2)
{
    return angle( v0-v1, v2-v1 );
}

/** Return the dihedral angle between v0-v1-v2-v3 (treating the vectors as points) */
MultiDouble MultiVector::dihedral(const MultiVector &v0, const MultiVector &v1,
                                  const MultiVector &v2, const MultiVector &v3)
{
    MultiDouble result;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        result.set( i, Vector::dihedral(v0.at(i),v1.at(i),v2.at(i),v3.at(i)) );
    }
    
    return result;
}
    
/** Generate a vector, v0, that has distance 'dst' v0-v1, angle 'ang' v0-v1-v2,
    and dihedral 'dih' v0-v1-v2-v3 */
MultiVector MultiVector::generate(
                   const MultiDouble &dst, const MultiVector &v1, const MultiDouble &ang,
                   const MultiVector &v2, const MultiDouble &dih, const MultiVector &v3)
{
    MultiVector result;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        result.set(i, Vector::generate(dst.at(i), v1.at(i), Angle(ang.at(i)),
                   v2.at(i), Angle(dih.at(i)), v3.at(i)));
    }
    
    return result;
}

/** Return the cross product of v0 and v1 */
MultiVector MultiVector::cross(const MultiVector &v0, const MultiVector &v1)
{
    MultiDouble nx = v0.sc[1]*v1.sc[2] - v0.sc[2]*v1.sc[1];
    MultiDouble ny = v0.sc[2]*v1.sc[0] - v0.sc[0]*v1.sc[2];
    MultiDouble nz = v0.sc[0]*v1.sc[1] - v0.sc[1]*v1.sc[0];

    MultiDouble length(nx);
    length *= nx;
    length.multiplyAdd(ny, ny);
    length.multiplyAdd(nz, nz);
    length = length.sqrt();
    
    MultiDouble near_parallel = length.compareLess( MultiDouble(0.01) );

    if (near_parallel.sum() > 0)
    {
        qDebug() << "NEAR PARALLEL" << near_parallel.toString()
                 << length.toString();
    
        //at least one of the pairs of vectors is parallel (or near parallel)
        // - manually calculate each cross product
        int nfixed = 0;
        
        for (int i=0; i<MultiDouble::count(); ++i)
        {
            if (near_parallel.at(i) > 0)
            {
                Vector normal = Vector::cross(v0.at(i), v1.at(i));
                nx.set(i, normal.x());
                ny.set(i, normal.y());
                nz.set(i, normal.z());
                length.set(i, 1.0);
                nfixed += 1;
            }
        }
        
        if (nfixed == MultiDouble::count())
        {
            //we have already made the vector
            return MultiVector(nx, ny, nz);
        }
    }

    //return the normalised vector
    MultiDouble mask = length.compareEqual( MultiDouble(0) );
    MultiDouble inv_length = mask.logicalAnd(length.reciprocal());
    
    return MultiVector( nx*inv_length, ny*inv_length, nz*inv_length );
}

/** Return the manhattan length of the vector */
MultiDouble MultiVector::manhattanLength() const
{
    MultiDouble result;
    
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        result.set(i, this->at(i).manhattanLength());
    }

    return result;
}

/** Increment, decrement, negate etc. */
MultiVector& MultiVector::operator/=(const MultiDouble &val)
{
    for (int i=0; i<MultiDouble::count(); ++i)
    {
        if (SireMaths::isZero(val.at(i)))
            throw SireMaths::math_error(QObject::tr(
                "Cannot divide a vector by zero! %1 / %2 is a error!")
                    .arg(this->toString()).arg(val.toString()),CODELOC);
    }

    for (int i=0; i<3; i++)
        sc[i] /= val;

    return *this;
}

namespace SireMaths
{
    /** Increment, decrement, negate etc. */
    MultiVector SIREMATHS_EXPORT operator/(const MultiVector &p1, const MultiDouble &c)
    {
        MultiVector result(p1);
        result /= c;
        return result;
    }
}

const char* MultiVector::typeName()
{
    return "SireMaths::MultiVector";
}

/** Swap the values of the value at index idx0 in 'f0' with the value at index 'idx' in 'f1' */
void MultiVector::swap(MultiVector &v0, int idx0, MultiVector &v1, int idx1)
{
    MultiDouble::swap(v0.sc[0], idx0, v1.sc[0], idx1);
    MultiDouble::swap(v0.sc[1], idx0, v1.sc[1], idx1);
    MultiDouble::swap(v0.sc[2], idx0, v1.sc[2], idx1);
}
