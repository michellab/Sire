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

#include <QDataStream>

#include "SireStream/datastream.h"

#include "axisset.h"

using namespace SireMaths;
using namespace SireStream;

static const RegisterMetaType<AxisSet> r_axisset(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const AxisSet &ax)
{
    writeHeader(ds, r_axisset, 1) << ax.mat << ax.invmat << ax.orgn;
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, AxisSet &ax)
{
    VersionID v = readHeader(ds, r_axisset);

    if (v == 1)
    {
        ds >> ax.mat >> ax.invmat >> ax.orgn;
    }
    else
        throw version_error(v, "1", r_axisset, CODELOC);

    return ds;
}

/** Construct an empty AxisSet. This represents the standard cartesian axes, centered
    on the origin */
AxisSet::AxisSet() : mat(), invmat(), orgn()
{}

/** Construct an AxisSet using matrix 'matrx', and origin 'orign' */
AxisSet::AxisSet(const Matrix &matrx, Vector vec)
        : mat(matrx), invmat(mat.inverse()), orgn(vec)
{}

/** Copy constructor */
AxisSet::AxisSet(const AxisSet &other)
        : mat(other.mat), invmat(other.invmat), orgn(other.orgn)
{}

/** Destructor */
AxisSet::~AxisSet()
{}

/** Convert a vector from the cartesian frame to this coordinate frame */
Vector AxisSet::fromIdentity(const Vector &vec) const
{
    return (mat * vec) + orgn;
}

/** Convert a vector from the cartesian frame with origin 'delta' to this coordinate frame */
Vector AxisSet::fromIdentity(const Vector &vec, const Vector &delta) const
{
    return delta + this->fromIdentity(vec-delta);
}

/** Convert the array of vectors from the cartesian frame to this coordinate frame */
QVector<Vector> AxisSet::fromIdentity(const QVector<Vector> &vecs) const
{
    QVector<Vector> newvecs(vecs);
    
    for (int i=0; i<vecs.count(); ++i)
    {
        newvecs[i] = this->fromIdentity(vecs[i]);
    }
    
    return newvecs;
}

/** Convert the array of vectors from the cartesian frame offset by delta
    to this coordinate frame */
QVector<Vector> AxisSet::fromIdentity(const QVector<Vector> &vecs, const Vector &delta) const
{
    QVector<Vector> newvecs(vecs);
    
    for (int i=0; i<vecs.count(); ++i)
    {
        newvecs[i] = this->fromIdentity(vecs[i],delta);
    }
    
    return newvecs;
}

/** Convert a vector to the cartesian frame from this coordinate frame */
Vector AxisSet::toIdentity(const Vector &vec) const
{
    return invmat * (vec - orgn);
}

/** Convert a vector from the frame 'frame' to this coordinate frame */
Vector AxisSet::fromFrame(const AxisSet &frame, const Vector &vec) const
{
    return ( mat *  (frame.invmat * (vec - frame.orgn)) ) + orgn;
}

/** Convert a vector to the frame 'frame' from this coordinate frame */
Vector AxisSet::toFrame(const AxisSet &frame, const Vector &vec) const
{
    return ( frame.mat * (invmat * (vec - orgn)) ) + frame.orgn;
}

/** Return a string representation of the AxisSet */
QString AxisSet::toString() const
{
    return QObject::tr("AxisSet{ origin() = %1,matrix=\n%2 }").arg(orgn.toString(),mat.toString());
}

const char* AxisSet::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AxisSet>() );
}
