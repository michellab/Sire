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

#ifndef SIREMOL_ATOMCOORDS_H
#define SIREMOL_ATOMCOORDS_H

#include "atomproperty.hpp"

#include "SireVol/coordgroup.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<>
class AtomProperty<SireMaths::Vector>;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::AtomProperty<SireMaths::Vector>&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::AtomProperty<SireMaths::Vector>&);

namespace SireMaths
{
class Vector;
class Quaternion;
class AxisSet;
class Matrix;
class Transform;
}

namespace SireVol
{
class Space;
}

namespace SireMol
{

using SireMaths::Vector;
using SireMaths::Quaternion;
using SireMaths::AxisSet;
using SireMaths::Matrix;
using SireMaths::Transform;

using SireVol::CoordGroup;
using SireVol::CoordGroupArray;
using SireVol::Space;

/** This is an explicit specialisation of AtomProperty<T> for the Vector
    class, as the Vector implies coordinates, which are arranged into
    CoordGroups (so that bounding boxes are calculated automatically)

    @author Christopher Woods
*/
template<>
class SIREMOL_EXPORT AtomProperty<Vector>
       : public SireBase::ConcreteProperty<AtomProperty<Vector>,AtomProp>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const AtomProperty<SireMaths::Vector>&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, AtomProperty<SireMaths::Vector>&);

public:
    AtomProperty();

    AtomProperty(const MoleculeInfoData &molinfo);

    AtomProperty(const CoordGroup &cgroup);
    AtomProperty(const CoordGroupArray &cgroups);

    AtomProperty(const AtomProperty<Vector> &other);

    ~AtomProperty();

    AtomProperty<Vector>& operator=(const AtomProperty<Vector> &other);

    static const char* typeName();

    AtomProperty<Vector>* clone() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    bool operator==(const AtomProperty<Vector> &other) const;
    bool operator!=(const AtomProperty<Vector> &other) const;

    bool canConvert(const QVariant &value) const;

    void assignFrom(const AtomProperty<QVariant> &values);

    AtomProperty<QVariant> toVariant() const;

    static AtomProperty<Vector> fromVariant(const AtomProperty<QVariant> &variant);

    const CoordGroup& operator[](CGIdx cgidx) const;

    const CoordGroup& at(CGIdx cgidx) const;
    const CoordGroup& get(CGIdx cgidx) const;

    const Vector& operator[](const CGAtomIdx &cgatomidx) const;
    const Vector& at(const CGAtomIdx &cgatomidx) const;
    const Vector& get(const CGAtomIdx &cgatomidx) const;

    AtomProperty<Vector>& set(const CGAtomIdx &cgatomidx, const Vector &value);

    AtomProperty<Vector>& set(CGIdx cgidx, const QVector<Vector> &values);
    AtomProperty<Vector>& set(CGIdx cgidx, const CoordGroup &cgroup);

    void translate(const Vector &delta);
    void translate(CGIdx cgidx, const Vector &delta);
    
    void rotate(const Quaternion &quat, const Vector &point);
    void rotate(const Matrix &rotmat, const Vector &point);
    
    void rotate(CGIdx cgidx, const Quaternion &quat, const Vector &point);
    void rotate(CGIdx cgidx, const Matrix &rotmat, const Vector &point);
    
    void transform(const Transform &t);
    void transform(CGIdx cgidx, const Transform &t);
    
    void mapInto(const AxisSet &axes);
    void mapInto(CGIdx cgidx, const AxisSet &axes);
    
    void changeFrame(const AxisSet &from_frame, const AxisSet &to_frame);
    void changeFrame(CGIdx cgidx, const AxisSet &from_frame, 
                                  const AxisSet &to_frame);

    const CoordGroup* data() const;
    const CoordGroup* constData() const;

    const Vector* data(CGIdx cgidx) const;
    const Vector* constData(CGIdx cgidx) const;

    int size() const;
    int count() const;

    int nCutGroups() const;

    int nAtoms() const;
    int nAtoms(CGIdx cgidx) const;

    const CoordGroupArray& array() const;

    QVector<Vector> toVector() const;
    QVector<Vector> toVector(const AtomSelection &selection) const;
    
    void copyFrom(const QVector<Vector> &values);
    void copyFrom(const QVector<Vector> &values, const AtomSelection &selection);

    PropertyPtr merge(const MoleculeInfoData &molinfo) const;
    PropertyPtr divide(const QVector<AtomSelection> &beads) const;
    PropertyPtr divideByResidue(const MoleculeInfoData &molinfo) const;

    void assertCanConvert(const QVariant &value) const;

private:
    /** The actual atomic coordinates, arranged into CoordGroups */
    CoordGroupArray coords;
};

/** Return the raw array that hold the coordinates */
inline const CoordGroupArray& AtomProperty<Vector>::array() const
{
    return coords;
}

typedef AtomProperty<Vector> AtomCoords;

}

Q_DECLARE_METATYPE( SireMol::AtomCoords );

SIRE_EXPOSE_ATOM_PROPERTY( SireMaths::Vector, SireMol::AtomCoords )

SIRE_END_HEADER

#endif
