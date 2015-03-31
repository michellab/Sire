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

#include "SireMaths/align.h"
#include "SireMaths/accumulator.h"

#include "SireStream/datastream.h"

#include "SireError/errors.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <QElapsedTimer>

#include "tostring.h"

using namespace SireMaths;
using namespace SireStream;

/////////
///////// Implementation of Transform
/////////

static const RegisterMetaType<Transform> r_trans(NO_ROOT);

QDataStream SIREMATHS_EXPORT &operator<<(QDataStream &ds, const Transform &trans)
{
    writeHeader(ds, r_trans, 1);
    
    ds << trans.delta << trans.rotcent << trans.rotmat;
    
    return ds;
}

QDataStream SIREMATHS_EXPORT &operator>>(QDataStream &ds, Transform &trans)
{
    VersionID v = readHeader(ds, r_trans);
    
    if (v == 1)
    {
        ds >> trans.delta >> trans.rotcent >> trans.rotmat;
    }
    else
        throw version_error(v, "1", r_trans, CODELOC);
    
    return ds;
}

/** Constructor (no transformation) */
Transform::Transform() : delta(0), rotcent(0), rotmat( Quaternion::identity() )
{}

/** Construct to only translate by 'delta' */
Transform::Transform(const Vector &del)
          : delta(del), rotcent(0), rotmat( Quaternion::identity() )
{}

/** Construct to only rotate using 'rotmat' around the rotation center 'center' */
Transform::Transform(const Quaternion &mat, const Vector &center)
          : delta(0), rotcent(center), rotmat(mat)
{
    if (rotmat.isIdentity())
    {
        rotcent = Vector(0);
    }
}

/** Construct to only rotate using 'rotmat' around the rotation center 'center' */
Transform::Transform(const Matrix &mat, const Vector &center)
          : delta(0), rotcent(center), rotmat(mat)
{
    if (rotmat.isIdentity())
    {
        rotcent = Vector(0);
    }
}

/** Construct to translate by delta and to rotate using 'rotmat' around the rotation
    center 'center' */
Transform::Transform(const Vector &del, const Quaternion &mat, const Vector &center)
          : delta(del), rotcent(center), rotmat(mat)
{
    if (rotmat.isIdentity())
    {
        rotcent = Vector(0);
    }
}

/** Construct to translate by delta and to rotate using 'rotmat' around the rotation
    center 'center' */
Transform::Transform(const Vector &del, const Matrix &mat, const Vector &center)
          : delta(del), rotcent(center), rotmat(mat)
{
    if (rotmat.isIdentity())
    {
        rotcent = Vector(0);
    }
}

/** Copy constructor */
Transform::Transform(const Transform &other)
          : delta(other.delta), rotcent(other.rotcent), rotmat(other.rotmat)
{}

/** Destructor */
Transform::~Transform()
{}

/** Copy assignment operator */
Transform& Transform::operator=(const Transform &other)
{
    if (this != &other)
    {
        delta = other.delta;
        rotcent = other.rotcent;
        rotmat = other.rotmat;
    }
    
    return *this;
}

/** Comparison operator */
bool Transform::operator==(const Transform &other) const
{
    return delta == other.delta and rotcent == other.rotcent and rotmat == other.rotmat;
}

/** Comparison operator */
bool Transform::operator!=(const Transform &other) const
{
    return not operator==(other);
}

const char* Transform::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Transform>() );
}

const char* Transform::what() const
{
    return Transform::typeName();
}

/** Return whether this is null (has no transformation) */
bool Transform::isNull() const
{
    return delta.isZero() and rotmat.isIdentity();
}

/** Return whether this is zero (has no transformation) */
bool Transform::isZero() const
{
    return delta.isZero() and rotmat.isIdentity();
}

QString Transform::toString() const
{
    if (isZero())
        return QObject::tr("Transform::null");
    
    else if (delta.isZero())
    {
        return QObject::tr("Transform( rotate by %1 about %2 )")
                    .arg(rotmat.toString()).arg(rotcent.toString());
    }
    else if (rotmat.isIdentity())
    {
        return QObject::tr("Transform( translate by %1 )")
                    .arg(delta.toString());
    }
    else
    {
        return QObject::tr("Transform( rotate by %1 about %2 then translate by %3 )")
                .arg(rotmat.toString()).arg(rotcent.toString())
                .arg(delta.toString());
    }
}

/** Apply this transformation to the passed point, returning the result */
Vector Transform::operator()(const Vector &point) const
{
    return delta + rotcent + rotmat.rotate(point-rotcent);
}

/** Apply this transformation to all of the passed points. Note that this
    modifies the points themselves and returns a pointer to the same array */
Vector* Transform::apply(Vector *points, int sz) const
{
    if (sz == 0 or (delta.isZero() and rotmat.isIdentity()))
        return points;
    
    else
    {
        if (rotmat.isIdentity())
        {
            for (int i=0; i<sz; ++i)
            {
                points[i] += delta;
            }
        }
        else if (delta.isZero())
        {
            if (rotcent.isZero())
            {
                for (int i=0; i<sz; ++i)
                {
                    points[i] = rotmat.rotate(points[i]);
                }
            }
            else
            {
                for (int i=0; i<sz; ++i)
                {
                    points[i] = rotcent + rotmat.rotate(points[i]-rotcent);
                }
            }
        }
        else if (rotcent.isZero())
        {
            for (int i=0; i<sz; ++i)
            {
                points[i] = delta + rotmat.rotate(points[i]);
            }
        }
        else
        {
            for (int i=0; i<sz; ++i)
            {
                points[i] = delta + rotcent + rotmat.rotate(points[i]-rotcent);
            }
        }
        
        return points;
    }
}

/** Apply this transformation to all of the passed points, returning the results */
QVector<Vector> Transform::operator()(const QVector<Vector> &points) const
{
    if (points.isEmpty() or (delta.isZero() and rotmat.isIdentity()))
        return points;
    else
    {
        QVector<Vector> ret(points);
        this->apply(ret.data(), ret.count());
        return ret;
    }
}

/** Apply this transformation to the passed point, returning the result */
Vector Transform::apply(const Vector &point) const
{
    return this->operator()(point);
}

/** Apply this transformation to all of the passed points, returning the results */
QVector<Vector> Transform::apply(const QVector<Vector> &points) const
{
    return this->operator()(points);
}

/** Return the amount by which to translate */
Vector Transform::translationDelta() const
{
    return delta;
}

/** Return the center of rotation */
Vector Transform::rotationCenter() const
{
    return rotcent;
}

/** Return the rotation matrix as a quaternion */
Quaternion Transform::rotationQuaternion() const
{
    return rotmat;
}

/** Return the rotation matrix */
Matrix Transform::rotationMatrix() const
{
    return rotmat.toMatrix();
}

/////////
///////// Alignment functions
/////////

static QVector<Vector> translate(const QVector<Vector> &v, const Vector &delta)
{
    QVector<Vector> v2(v);
    
    for (int i=0; i<v.count(); ++i)
    {
        v2[i] += delta;
    }
    
    return v2;
}

static QVector<Vector> rotate(const QVector<Vector> &v, const Matrix &rotmat)
{
    QVector<Vector> v2(v);
    
    for (int i=0; i<v.count(); ++i)
    {
        v2[i] = rotmat * v2[i];
    }
    
    return v2;
}

namespace SireMaths
{
    /** Return the centroid of the points in 'p'. If n != -1 then
        only calculate the centroid of the first n points */
    Vector SIREMATHS_EXPORT getCentroid(const QVector<Vector> &p, int n)
    {
        if (p.isEmpty())
            return Vector(0);
    
        Average x, y, z;
    
        if (n < 0 or n > p.count())
        {
            n = p.count();
        }
        
        for (int i=0; i<n; ++i)
        {
            x.accumulate(p[i].x());
            y.accumulate(p[i].y());
            z.accumulate(p[i].z());
        }
        
        return Vector(x.average(), y.average(), z.average());
    }

    /** Return the RMSD between the two sets of points. If n != -1 then
        only calculate the RMSD using the first n points */
    double SIREMATHS_EXPORT getRMSD(const QVector<Vector> &p, const QVector<Vector> &q, int n)
    {
        if (p.isEmpty() or q.isEmpty())
            return 0;
        
        if (n < 0 or n > qMin(p.count(),q.count()))
        {
            n = qMin(p.count(), q.count());
        }
        
        Average msd;
        
        for (int i=0; i<n; ++i)
        {
            msd.accumulate( pow_2(p[i].x()-q[i].x()) +
                            pow_2(p[i].y()-q[i].y()) +
                            pow_2(p[i].z()-q[i].z()) );
        }
        
        return std::sqrt( msd.average() );
    }

    /** Use the kabasch algorithm (http://en.wikipedia.org/wiki/Kabsch_algorithm)
        to calculate the rotation matrix needed to align
        'q' on top of 'p'. This will match p[i] against
        q[i], and only matches the first N points (where N is the minimum
        of the number of points in 'p' and the number of points in 'q'.
        
        Note that this code assumes that q has been already translated on top of p.
     
        This code is inspired by the calculate_rmsd python script written by
        Jimmy Charnley Kromann and Lars Bratholm, available
        https://github.com/charnley/rmsd, and under license;
        
        =====================
        Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice, this
           list of conditions and the following disclaimer.
        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        =======================
     
        (note that the C++ implement that I have here is licensed under the GPL,
         as stated at the top of this file)
    */
    Matrix SIREMATHS_EXPORT kabasch(const QVector<Vector> &p,
                                    const QVector<Vector> &q)
    {
        if (p.isEmpty() or q.isEmpty())
            return Matrix::identity();

        //calculate the covariance matrix
        Matrix c = Matrix::covariance(p, q);
        
        //calculate the single value decomposition of this in V S W^T
        boost::tuple<Matrix,Matrix,Matrix> svd = c.svd();
        
        Matrix v = svd.get<0>();
        Matrix s = svd.get<1>();
        Matrix w = svd.get<2>();
        
        double det_vw = v.determinant() * w.determinant();

        if (det_vw < 0)
        {
            for (int i=0; i<3; ++i)
            {
                v(i,2) = -v(i,2);
            }
        }

        //now create the rotation matrix
        Matrix r = v * w;

        //ensure that this is a rotation matrix (has a determinant of 1)
        const double det = r.determinant();
        
        if (not SireMaths::areEqual(det, 1.0))
        {
            throw SireError::program_bug( QObject::tr(
                    "Attempt to find the alignment matrix failed to produce a valid "
                    "rotation matrix. Determinant should be 1. It is equal to %1.")
                        .arg(det), CODELOC );
        }

        return r;
    }
    
    /** Use the kabasch algorithm (http://en.wikipedia.org/wiki/Kabsch_algorithm)
        to calculate the translation vector and rotation matrix needed to align
        'q' on top of 'p'. This will match p[i] against
        q[i], and only matches the first N points (where N is the minimum
        of the number of points in 'p' and the number of points in 'q'.
        This will provide additional fitting to find the optimum translation vector
        for the alignment (as opposed to merely using the difference between the 
        centroids, as used in the simple kabasch algorithm)
     
        This code is inspired by the calculate_rmsd python script written by
        Jimmy Charnley Kromann and Lars Bratholm, available
        https://github.com/charnley/rmsd, and under license;
        
        =====================
        Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

        1. Redistributions of source code must retain the above copyright notice, this
           list of conditions and the following disclaimer.
        2. Redistributions in binary form must reproduce the above copyright notice,
           this list of conditions and the following disclaimer in the documentation
           and/or other materials provided with the distribution.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
        ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
        WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
        ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
        (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
        LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
        (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        =======================
     
        (note that the C++ implement that I have here is licensed under the GPL,
         as stated at the top of this file)
    */
    Transform SIREMATHS_EXPORT kabaschFit(const QVector<Vector> &p,
                                          const QVector<Vector> &q)
    {
        if (p.isEmpty() or q.isEmpty())
            return Transform();
        
        //find the maximum point in p
        Vector step_size(0);
        
        for (int i=0; i<qMin(p.count(),q.count()); ++i)
        {
            step_size = step_size.max(q[i]);
        }
        
        const Vector threshold = 1e-9 * step_size;
        
        QVector<Vector> q2(q);
    
        q2 = translate(q, Vector(0,0,0) );
    
        Matrix rotmat_best = kabasch(p,q2);
        
        double rmsd_best = getRMSD( p, rotate(q2,rotmat_best) );
        
        Vector delta(0);
        
        while (true)
        {
            for (int i=0; i<3; ++i)
            {
                Vector tmp = delta;
                tmp.set(i, tmp[i] + step_size[i]);
                
                q2 = translate(q,tmp);
                Matrix rotmat = kabasch(p,q2);
                double rmsd = getRMSD( p, rotate(q2,rotmat) );
                
                if (rmsd < rmsd_best)
                {
                    //this has improved the fit
                    rmsd_best = rmsd;
                    rotmat_best = rotmat;
                    delta = tmp;
                }
                else
                {
                    //try the other way
                    tmp = delta;
                    tmp.set(i, tmp[i] - step_size[i]);
                    
                    q2 = translate(q,tmp);
                    rotmat = kabasch(p,q2);
                    rmsd = getRMSD( p, rotate(q2,rotmat) );
                    
                    if (rmsd < rmsd_best)
                    {
                        //this has improved the fit
                        rmsd_best = rmsd;
                        rotmat_best = rotmat;
                        delta = tmp;
                    }
                    else
                    {
                        //we need to reduce the amount by which we are searching
                        step_size.set(i, 0.5*step_size[i]);
                    }
                }
            }
            
            bool finished = true;
            
            for (int i=0; i<3; ++i)
            {
                if (step_size[i] > threshold[i])
                {
                    finished = false;
                    break;
                }
            }
            
            if (finished)
                break;
        }
        
        return Transform(delta, rotmat_best, Vector(0));
    }

    /** Return AxisSet containing the translation vector and rotation
        matrix needed to optimally align the points in 'q' on top of the
        points in 'p'. If 'fit' is true, then this performs an RMSD
        fit to find the optimal translation vector (as opposed to merely
        taking the difference of centroids) */
    Transform SIREMATHS_EXPORT getAlignment(const QVector<Vector> &p,
                                            const QVector<Vector> &q,
                                            bool fit)
    {
        if (p.isEmpty() or q.isEmpty())
            return Transform();

        //first, translate p and q so that their centroids
        //are at (0,0,0)

        //calculate the difference in centroids of p and q
        const int n = qMin(p.count(), q.count());
        
        Vector cp = getCentroid(p,n);
        Vector cq = getCentroid(q,n);
        
        QVector<Vector> pc(n), qc(n);
        for (int i=0; i<n; ++i)
        {
            pc[i] = p[i] - cp;
            qc[i] = q[i] - cq;
        }
    
        if (fit)
        {
            Transform a = kabaschFit(pc, qc);
            
            return Transform(a.translationDelta() + cp - cq, a.rotationQuaternion(), cq);
        }
        else
        {
            Matrix rotmat = kabasch(pc, qc);
            
            return Transform(cp-cq, rotmat, cq);
        }
    }
    
    /** Return a copy of points 'q' aligned on top of points 'p'. If 'fit' 
        is true, then this performs an RMSD
        fit to find the optimal translation vector (as opposed to merely
        taking the difference of centroids) */
    QVector<Vector> SIREMATHS_EXPORT align(const QVector<Vector> &p,
                                           const QVector<Vector> &q,
                                           bool fit)
    {
        if (p.isEmpty() or q.isEmpty())
            return q;
        
        Transform a = getAlignment(p, q, fit);
        
        return a.apply(q);
    }

} // end of namespace SireMaths
