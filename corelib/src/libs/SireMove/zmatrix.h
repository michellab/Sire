/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#ifndef SIREMOVE_ZMATRIX_H
#define SIREMOVE_ZMATRIX_H

#include "SireBase/propertymap.h"

#include "SireMol/molviewproperty.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomidx.h"

#include <QHash>
#include <QVector>

SIRE_BEGIN_HEADER

namespace SireMove
{
class ZMatrix;
class ZMatrixLine;

class ZMatrixCoords;
class ZMatrixCoordsLine;
}

QDataStream& operator<<(QDataStream &ds, const SireMove::ZMatrix&);
QDataStream& operator>>(QDataStream &ds, SireMove::ZMatrix&);

QDataStream& operator<<(QDataStream &ds, const SireMove::ZMatrixLine&);
QDataStream& operator>>(QDataStream &ds, SireMove::ZMatrixLine&);

QDataStream& operator<<(QDataStream &ds, const SireMove::ZMatrixCoords&);
QDataStream& operator>>(QDataStream &ds, SireMove::ZMatrixCoords&);

QDataStream& operator<<(QDataStream &ds, const SireMove::ZMatrixCoordsLine&);
QDataStream& operator>>(QDataStream &ds, SireMove::ZMatrixCoordsLine&);

namespace SireMol
{
class PartialMolecule;
class BondID;
class AngleID;
class DihedralID;
}

namespace SireMove
{

using SireBase::PropertyMap;

using SireUnits::Dimension::Length;
using SireUnits::Dimension::Angle;

using SireMol::AtomCoords;
using SireMol::AtomIdx;
using SireMol::AtomID;
using SireMol::BondID;
using SireMol::AngleID;
using SireMol::DihedralID;
using SireMol::PartialMolecule;
using SireMol::Molecule;
using SireMol::MoleculeInfoData;
using SireMol::AtomSelection;
using SireMol::AtomMatcher;

using SireMaths::Vector;

/** This class represents a single line of a zmatrix

    @author Christopher Woods
*/
class SIREMOVE_EXPORT ZMatrixLine
{

friend QDataStream& ::operator<<(QDataStream&, const ZMatrixLine&);
friend QDataStream& ::operator>>(QDataStream&, ZMatrixLine&);

public:
    ZMatrixLine();
    ZMatrixLine(AtomIdx atom, AtomIdx bond,
                AtomIdx angle, AtomIdx dihedral);

    ZMatrixLine(const ZMatrixLine &other);

    ~ZMatrixLine();

    ZMatrixLine& operator=(const ZMatrixLine &other);

    bool operator==(const ZMatrixLine &other) const;
    bool operator!=(const ZMatrixLine &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ZMatrixLine::typeName();
    }

    QString toString() const;

    AtomIdx operator[](int i) const;

    AtomIdx atom() const;
    AtomIdx bond() const;
    AtomIdx angle() const;
    AtomIdx dihedral() const;

    Length bondDelta() const;
    Angle angleDelta() const;
    Angle dihedralDelta() const;

    void setBondDelta(const Length &delta);
    void setAngleDelta(const Angle &delta);
    void setDihedralDelta(const Angle &delta);

private:
    /** The indicies of the four atoms */
    AtomIdx atms[4];

    /** The bond, angle and dihedral deltas */
    Vector deltas;
};

/** This class holds a z-matrix line that includes the
    sizes of the internal coordinates

    @author Christopher Woods
*/
class SIREMOVE_EXPORT ZMatrixCoordsLine : public ZMatrixLine
{

friend QDataStream& ::operator<<(QDataStream&, const ZMatrixCoordsLine&);
friend QDataStream& ::operator>>(QDataStream&, ZMatrixCoordsLine&);

public:
    ZMatrixCoordsLine();
    ZMatrixCoordsLine(const ZMatrixLine &line);
    ZMatrixCoordsLine(const ZMatrixLine &line,
                      const Length &bond, const Angle &angle, const Angle &dihedral);

    ZMatrixCoordsLine(const ZMatrixCoordsLine &other);

    ~ZMatrixCoordsLine();

    ZMatrixCoordsLine& operator=(const ZMatrixCoordsLine &other);

    bool operator==(const ZMatrixCoordsLine &other) const;
    bool operator!=(const ZMatrixCoordsLine &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ZMatrixCoordsLine::typeName();
    }

    QString toString() const;

    Length bondLength() const;
    Angle angleSize() const;
    Angle dihedralSize() const;

    void setBond(const Length &length);
    void setAngle(const Angle &size);
    void setDihedral(const Angle &size);

private:
    /** The internal coordinates */
    Vector coords;
};

/** This class holds the z-matrix of a molecule. A z-matrix
    is be used to provide an internal coordinates representation
    of a molecule, so that internal moves may be performed

    @author Christopher Woods
*/
class SIREMOVE_EXPORT ZMatrix
          : public SireBase::ConcreteProperty<ZMatrix,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const ZMatrix&);
friend QDataStream& ::operator>>(QDataStream&, ZMatrix&);

public:
    ZMatrix();

    ZMatrix(const Molecule &molecule);
    ZMatrix(const MoleculeInfoData &molinfo);

    ZMatrix(const ZMatrix &other);

    ~ZMatrix();

    ZMatrix& operator=(const ZMatrix &other);

    bool operator==(const ZMatrix &other) const;
    bool operator!=(const ZMatrix &other) const;

    static const char* typeName();

    const ZMatrixLine& operator[](const AtomID &atom) const;

    const ZMatrixLine& at(const AtomID &atom) const;

    const QVector<ZMatrixLine>& lines() const;

    const QHash<AtomIdx,int>& index() const;

    int nLines() const;

    int getIndex(const AtomID &atom) const;
    int getIndex(const AtomID &atom, const AtomID &bond) const;
    int getIndex(const AtomID &atom, const AtomID &bond,
                 const AtomID &angle) const;
    int getIndex(const AtomID &atom, const AtomID &bond,
                 const AtomID &angle, const AtomID &dihedral) const;

    int getIndex(const BondID &bond) const;
    int getIndex(const AngleID &angle) const;
    int getIndex(const DihedralID &dihedral) const;

    QString toString() const;

    const SireMol::MoleculeInfoData& info() const;

    bool contains(const AtomID &atom) const;
    bool contains(const AtomID &atom, const AtomID &bond) const;
    bool contains(const AtomID &atom, const AtomID &bond,
                  const AtomID &angle) const;
    bool contains(const AtomID &atom, const AtomID &bond,
                  const AtomID &angle, const AtomID &dihedral) const;

    bool contains(const BondID &bond) const;
    bool contains(const AngleID &angle) const;
    bool contains(const DihedralID &dihedral) const;

    void add(const AtomID &atom, const AtomID &bond,
             const AtomID &angle, const AtomID &dihedral);

    void add(const DihedralID &dihedral);

    void remove(const AtomID &atom);
    void remove(const AtomID &atom, const AtomID &bond,
                const AtomID &angle, const AtomID &dihedral);

    void remove(const DihedralID &dihedral);

    void add(const ZMatrixLine &zmatline);
    void remove(const ZMatrixLine &zmatline);

    void setDelta(const BondID &bond, const Length &delta);
    void setDelta(const AngleID &angle, const Angle &delta);
    void setDelta(const DihedralID &dihedral, const Angle &delta);

    void setBondDelta(const AtomID &atom, const Length &delta);
    void setAngleDelta(const AtomID &atom, const Angle &delta);
    void setDihedralDelta(const AtomID &atom, const Angle &delta);

    void setBondDelta(const AtomID &atom, const AtomID &bond,
                      const Length &delta);

    void setAngleDelta(const AtomID &atom, const AtomID &bond,
                       const AtomID &angle, const Angle &delta);

    void setDihedralDelta(const AtomID &atom, const AtomID &bond,
                          const AtomID &angle, const AtomID &dihedral,
                          const Angle &delta);

    Length delta(const BondID &bond) const;
    Angle delta(const AngleID &angle) const;
    Angle delta(const DihedralID &dihedral) const;

    Length bondDelta(const AtomID &atom) const;
    Angle angleDelta(const AtomID &atom) const;
    Angle dihedralDelta(const AtomID &atom) const;

    Length bondDelta(const AtomID &atom, const AtomID &bond) const;
    Angle angleDelta(const AtomID &atom, const AtomID &bond,
                     const AtomID &angle) const;
    Angle dihedralDelta(const AtomID &atom, const AtomID &bond,
                        const AtomID &angle, const AtomID &dihedral) const;

    ZMatrix matchToSelection(const AtomSelection &selection) const;

    const QVector<int>& atomBuildOrder() const;

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

protected:
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const AtomMatcher &atommatcher) const;
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const QHash<AtomIdx,AtomIdx> &map) const;

private:
    void rebuildOrder();
    void reindex();

    /** The layout of the molecule whose coordinates
        are represented in this zmatrix */
    SireBase::SharedDataPointer<SireMol::MoleculeInfoData> molinfo;

    /** The z-matrix, held in construction order */
    QVector<ZMatrixLine> zmat;

    /** The index of each atom in the zmatrix */
    QHash<AtomIdx,int> atomidx_to_zmat;

    /** The order in which atoms should be constructed using
        this z-matrix */
    QVector<int> zmat_build_order;
};

/** This class holds a z-matrix of a molecule, together with the
    internal coordinate values for this z-matrix

    @author Christopher Woods
*/
class SIREMOVE_EXPORT ZMatrixCoords
            : public SireBase::ConcreteProperty<ZMatrixCoords,SireMol::MoleculeProperty>
{

friend QDataStream& ::operator<<(QDataStream&, const ZMatrixCoords&);
friend QDataStream& ::operator>>(QDataStream&, ZMatrixCoords&);

public:
    ZMatrixCoords();
    ZMatrixCoords(const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());

    ZMatrixCoords(const ZMatrix &zmatrix, const PartialMolecule &molecule,
                  const PropertyMap &map = PropertyMap());

    ZMatrixCoords(const ZMatrix &zmatrix, const AtomCoords &coords);

    ZMatrixCoords(const ZMatrixCoords &other);

    ~ZMatrixCoords();

    ZMatrixCoords& operator=(const ZMatrixCoords &other);

    bool operator==(const ZMatrixCoords &other) const;
    bool operator!=(const ZMatrixCoords &other) const;

    static const char* typeName();

    ZMatrixCoordsLine operator[](const AtomID &atom) const;

    ZMatrixCoordsLine at(const AtomID &atom) const;

    QVector<ZMatrixCoordsLine> lines() const;

    const QHash<AtomIdx,int>& index() const;

    int getIndex(const AtomID &atom) const;
    int getIndex(const AtomID &atom, const AtomID &bond) const;
    int getIndex(const AtomID &atom, const AtomID &bond,
                 const AtomID &angle) const;
    int getIndex(const AtomID &atom, const AtomID &bond,
                 const AtomID &angle, const AtomID &dihedral) const;

    int getIndex(const BondID &bond) const;
    int getIndex(const AngleID &angle) const;
    int getIndex(const DihedralID &dihedral) const;

    QString toString() const;

    const SireMol::MoleculeInfoData& info() const;

    const ZMatrix& zmatrix() const;

    bool contains(const AtomID &atom) const;

    bool contains(const AtomID &atom, const AtomID &bond) const;
    bool contains(const AtomID &atom, const AtomID &bond,
                  const AtomID &angle) const;
    bool contains(const AtomID &atom, const AtomID &bond,
                  const AtomID &angle, const AtomID &dihedral) const;

    bool contains(const BondID &bond) const;
    bool contains(const AngleID &angle) const;
    bool contains(const DihedralID &dihedral) const;

    void add(const AtomID &atom, const AtomID &bond,
             const AtomID &angle, const AtomID &dihedral);

    void add(const DihedralID &dihedral);

    void add(const AtomID &atom,
             const Length &bondlength, const AtomID &bond,
             const Angle &anglesize, const AtomID &angle,
             const Angle &dihedralsize, const AtomID &dihedral);

    void remove(const AtomID &atom);
    void remove(const AtomID &atom, const AtomID &bond,
                const AtomID &angle, const AtomID &dihedral);

    void remove(const DihedralID &dihedral);

    void add(const ZMatrixLine &zmatline);
    void remove(const ZMatrixLine &zmatline);

    void add(const ZMatrixCoordsLine &zmatline);

    bool isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const;

    const AtomCoords& toCartesian() const;

    const QVector<Vector>& internalCoordinates() const;

    void move(const BondID &bond, const Length &delta);
    void move(const AngleID &angle, const Angle &delta);
    void move(const DihedralID &dihedral, const Angle &delta);

    void moveBond(const AtomID &atom, const Length &delta);
    void moveAngle(const AtomID &atom, const Angle &delta);
    void moveDihedral(const AtomID &atom, const Angle &delta);

    void moveBond(const AtomID &atom0, const AtomID &atom1,
                  const Length &delta);

    void moveAngle(const AtomID &atom0, const AtomID &atom1,
                   const AtomID &atom2, const Angle &delta);

    void moveDihedral(const AtomID &atom0, const AtomID &atom1,
                      const AtomID &atom2, const AtomID &atom3,
                      const Angle &delta);

    void set(const BondID &bond, const Length &length);
    void set(const AngleID &angle, const Angle &size);
    void set(const DihedralID &dihedral, const Angle &size);

    void setBond(const AtomID &atom, const Length &length);
    void setAngle(const AtomID &atom, const Angle &size);
    void setDihedral(const AtomID &atom, const Angle &size);

    void setBond(const AtomID &atom0, const AtomID &atom1,
                 const Length &length);

    void setAngle(const AtomID &atom0, const AtomID &atom1,
                  const AtomID &atom2, const Angle &size);

    void setDihedral(const AtomID &atom0, const AtomID &atom1,
                     const AtomID &atom2, const AtomID &atom3,
                     const Angle &size);

    Length bondLength(const AtomID &atom) const;
    Angle angleSize(const AtomID &atom) const;
    Angle dihedralSize(const AtomID &atom) const;

    Length length(const BondID &bond) const;
    Angle size(const AngleID &angle) const;
    Angle size(const DihedralID &dihedral) const;

    Length bondLength(const AtomID &atom, const AtomID &bond) const;
    Angle angleSize(const AtomID &atom, const AtomID &bond,
                    const AtomID &angle) const;
    Angle dihedralSize(const AtomID &atom, const AtomID &bond,
                       const AtomID &angle, const AtomID &dihedral) const;

    void setDelta(const BondID &bond, const Length &delta);
    void setDelta(const AngleID &angle, const Angle &delta);
    void setDelta(const DihedralID &dihedral, const Angle &delta);

    void setBondDelta(const AtomID &atom, const Length &delta);
    void setAngleDelta(const AtomID &atom, const Angle &delta);
    void setDihedralDelta(const AtomID &atom, const Angle &delta);

    void setBondDelta(const AtomID &atom, const AtomID &bond,
                      const Length &delta);

    void setAngleDelta(const AtomID &atom, const AtomID &bond,
                       const AtomID &angle, const Angle &delta);

    void setDihedralDelta(const AtomID &atom, const AtomID &bond,
                          const AtomID &angle, const AtomID &dihedral,
                          const Angle &delta);

    Length delta(const BondID &bond) const;
    Angle delta(const AngleID &angle) const;
    Angle delta(const DihedralID &angle) const;

    Length bondDelta(const AtomID &atom) const;
    Angle angleDelta(const AtomID &atom) const;
    Angle dihedralDelta(const AtomID &atom) const;

    Length bondDelta(const AtomID &atom, const AtomID &bond) const;
    Angle angleDelta(const AtomID &atom, const AtomID &bond,
                     const AtomID &angle) const;
    Angle dihedralDelta(const AtomID &atom, const AtomID &bond,
                        const AtomID &angle, const AtomID &dihedral) const;

    ZMatrixCoords matchToSelection(const AtomSelection &selection) const;

private:
    void rebuildInternals();
    void rebuildCartesian() const;

    void _pvt_rebuildCartesian();

    void addInternal(const AtomIdx &atom);

    Vector getInternalCoords(const ZMatrixLine &line) const;

    /** The underlying z-matrix */
    ZMatrix zmat;

    /** The values of the internal coordinates */
    QVector<Vector> internal_coords;

    /** The current cartesian coordinates - these are used
        for the atoms that are not explicitly in the z-matrix */
    AtomCoords cartesian_coords;

    /** Whether or not the cartesian coordinates need to be rebuilt */
    bool need_rebuild;
};

}

Q_DECLARE_METATYPE( SireMove::ZMatrix )
Q_DECLARE_METATYPE( SireMove::ZMatrixLine )

Q_DECLARE_METATYPE( SireMove::ZMatrixCoords )
Q_DECLARE_METATYPE( SireMove::ZMatrixCoordsLine )

SIRE_EXPOSE_CLASS( SireMove::ZMatrix )
SIRE_EXPOSE_CLASS( SireMove::ZMatrixLine )

SIRE_EXPOSE_CLASS( SireMove::ZMatrixCoords )
SIRE_EXPOSE_CLASS( SireMove::ZMatrixCoordsLine )

SIRE_END_HEADER

#endif
