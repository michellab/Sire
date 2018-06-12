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

#include <QTime>

#include "zmatrix.h"

#include "SireID/index.h"

#include "SireMol/molecule.h"
#include "SireMol/partialmolecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/atommatcher.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"

#include "SireUnits/convert.h"
#include "SireUnits/units.h"

#include "SireMove/errors.h"
#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QTime>
#include <QDebug>

using namespace SireMove;
using namespace SireMol;
using namespace SireID;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of ZMatrixLine
//////////

static const RegisterMetaType<ZMatrixLine> r_zmatline(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const ZMatrixLine &zmatline)
{
    writeHeader(ds, r_zmatline, 1);

    for (int i=0; i<4; ++i)
        ds << zmatline.atms[i];

    ds << zmatline.deltas;

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, ZMatrixLine &zmatline)
{
    VersionID v = readHeader(ds, r_zmatline);

    if (v == 1)
    {
        for (int i=0; i<4; ++i)
            ds >> zmatline.atms[i];

        ds >> zmatline.deltas;
    }
    else
        throw version_error( v, "1", r_zmatline, CODELOC );

    return ds;
}

/** Null constructor */
ZMatrixLine::ZMatrixLine()
{
    for (int i=0; i<4; ++i)
        atms[i] = AtomIdx(0);
}

/** Construct with the specified atom indicies */
ZMatrixLine::ZMatrixLine(AtomIdx atom, AtomIdx bond,
                         AtomIdx angle, AtomIdx dihedral)
{
    atms[0] = atom;
    atms[1] = bond;
    atms[2] = angle;
    atms[3] = dihedral;
}

/** Copy constructor */
ZMatrixLine::ZMatrixLine(const ZMatrixLine &other)
            : deltas(other.deltas)
{
    for (int i=0; i<4; ++i)
        atms[i] = other.atms[i];
}

/** Destructor */
ZMatrixLine::~ZMatrixLine()
{}

/** Copy assignment operator */
ZMatrixLine& ZMatrixLine::operator=(const ZMatrixLine &other)
{
    if (this != &other)
    {
        deltas = other.deltas;

        for (int i=0; i<4; ++i)
            atms[i] = other.atms[i];
    }

    return *this;
}

/** Comparison operator */
bool ZMatrixLine::operator==(const ZMatrixLine &other) const
{
    return atms[0] == other.atms[0] and
           atms[1] == other.atms[1] and
           atms[2] == other.atms[2] and
           atms[3] == other.atms[3] and
           deltas == other.deltas;
}

/** Comparison operator */
bool ZMatrixLine::operator!=(const ZMatrixLine &other) const
{
    return not this->operator==(other);
}

/** Return a string representation */
QString ZMatrixLine::toString() const
{
    return QObject::tr("ZMatrix: %1-%2-%3-%4")
                 .arg(atms[0].value()).arg(atms[1].value())
                 .arg(atms[2].value()).arg(atms[3].value());
}

/** Return the ith atom index (i==0 is atom, i==1 is bond etc.)

    \throw SireError::invalid_index
*/
AtomIdx ZMatrixLine::operator[](int i) const
{
    return atms[ Index(i).map(4) ];
}

/** Return the index of the atom whose coordinates
    are described in this line */
AtomIdx ZMatrixLine::atom() const
{
    return atms[0];
}

/** Return the index of the bonded atom */
AtomIdx ZMatrixLine::bond() const
{
    return atms[1];
}

/** Return the index of the angled atom */
AtomIdx ZMatrixLine::angle() const
{
    return atms[2];
}

/** Return the index of the dihedralled atom */
AtomIdx ZMatrixLine::dihedral() const
{
    return atms[3];
}

/** Return the maximum amount by which the bond should be changed */
Length ZMatrixLine::bondDelta() const
{
    return Length( deltas[0] );
}

/** Return the maximum amount by which the angle should be changed */
Angle ZMatrixLine::angleDelta() const
{
    return Angle( deltas[1] );
}

/** Return the maximum amount by which the dihedral should be changed */
Angle ZMatrixLine::dihedralDelta() const
{
    return Angle( deltas[2] );
}

/** Set the maximum amount by which the bond should be changed */
void ZMatrixLine::setBondDelta(const Length &delta)
{
    deltas.setX( delta.value() );
}

/** Set the maximum amount by which the angle should be changed */
void ZMatrixLine::setAngleDelta(const Angle &delta)
{
    deltas.setY( delta.value() );
}

/** Set the maximum amount by which the dihedral should be changed */
void ZMatrixLine::setDihedralDelta(const Angle &delta)
{
    deltas.setZ( delta.value() );
}

const char* ZMatrixLine::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ZMatrixLine>() );
}

//////////
////////// Implementation of ZMatrixCoordsLine
//////////

static const RegisterMetaType<ZMatrixCoordsLine> r_zmatcoordsline(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds,
                                        const ZMatrixCoordsLine &zmatcoordsline)
{
    writeHeader(ds, r_zmatcoordsline, 1);

    ds << zmatcoordsline.coords
       << static_cast<const ZMatrixLine&>(zmatcoordsline);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds,
                                        ZMatrixCoordsLine &zmatcoordsline)
{
    VersionID v = readHeader(ds, r_zmatcoordsline);

    if (v == 1)
    {
        ds >> zmatcoordsline.coords
           >> static_cast<ZMatrixLine&>(zmatcoordsline);
    }
    else
        throw version_error( v, "1", r_zmatcoordsline, CODELOC );

    return ds;
}

/** Null constructor */
ZMatrixCoordsLine::ZMatrixCoordsLine() : ZMatrixLine()
{}

/** Construct from the passed line (but with zero coordinates) */
ZMatrixCoordsLine::ZMatrixCoordsLine(const ZMatrixLine &line)
                  : ZMatrixLine(line)
{}

/** Construct from the passed line and coordinates */
ZMatrixCoordsLine::ZMatrixCoordsLine(const ZMatrixLine &line,
                                     const Length &bond, const Angle &angle,
                                     const Angle &dihedral)
                  : ZMatrixLine(line)
{
    coords.setX( bond.value() );
    coords.setY( angle.value() );
    coords.setZ( dihedral.value() );
}

/** Copy constructor */
ZMatrixCoordsLine::ZMatrixCoordsLine(const ZMatrixCoordsLine &other)
                  : ZMatrixLine(other), coords(other.coords)
{}

/** Destructor */
ZMatrixCoordsLine::~ZMatrixCoordsLine()
{}

/** Copy assignment operator */
ZMatrixCoordsLine& ZMatrixCoordsLine::operator=(const ZMatrixCoordsLine &other)
{
    ZMatrixLine::operator=(other);
    coords = other.coords;

    return *this;
}

/** Comparison operator */
bool ZMatrixCoordsLine::operator==(const ZMatrixCoordsLine &other) const
{
    return ZMatrixLine::operator==(other) and coords == other.coords;
}

/** Comparison operator */
bool ZMatrixCoordsLine::operator!=(const ZMatrixCoordsLine &other) const
{
    return ZMatrixLine::operator!=(other) or coords != other.coords;
}

/** Return a string representation */
QString ZMatrixCoordsLine::toString() const
{
    return QObject::tr("%1 - %2 A : %3' : %4'")
                 .arg(ZMatrixLine::toString())
                 .arg(coords[0])
                 .arg( Angle(coords[1]).to(degrees) )
                 .arg( Angle(coords[2]).to(degrees) );
}

/** Return the length of the bond */
Length ZMatrixCoordsLine::bondLength() const
{
    return Length(coords[0]);
}

/** Return the size of the angle */
Angle ZMatrixCoordsLine::angleSize() const
{
    return Angle(coords[1]);
}

/** Return the size of the dihedral */
Angle ZMatrixCoordsLine::dihedralSize() const
{
    return Angle(coords[2]);
}

/** Set the length of the bond */
void ZMatrixCoordsLine::setBond(const Length &length)
{
    coords.setX( length.value() );
}

/** Set the size of the angle */
void ZMatrixCoordsLine::setAngle(const Angle &size)
{
    coords.setY( size.value() );
}

/** Set the size of the dihedral */
void ZMatrixCoordsLine::setDihedral(const Angle &size)
{
    coords.setZ( size.value() );
}

const char* ZMatrixCoordsLine::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ZMatrixCoordsLine>() );
}

//////////
////////// Implementation of ZMatrix
//////////

static const RegisterMetaType<ZMatrix> r_zmat;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const ZMatrix &zmat)
{
    writeHeader(ds, r_zmat, 1);

    SharedDataStream sds(ds);

    sds << zmat.molinfo << zmat.zmat
        << static_cast<const MoleculeProperty&>(zmat);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, ZMatrix &zmat)
{
    VersionID v = readHeader(ds, r_zmat);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> zmat.molinfo >> zmat.zmat
            >> static_cast<MoleculeProperty&>(zmat);

        zmat.reindex();
    }
    else
        throw version_error( v, "1", r_zmat, CODELOC );

    return ds;
}

/** Null constructor */
ZMatrix::ZMatrix() : ConcreteProperty<ZMatrix,MoleculeProperty>(),
                     molinfo( MoleculeInfoData::null() )
{}

/** Construct to hold the z-matrix for the passed molecule */
ZMatrix::ZMatrix(const Molecule &molecule)
        : ConcreteProperty<ZMatrix,MoleculeProperty>(),
          molinfo(molecule.data().info())
{}

/** Construct to hold the z-matrix for the molecule whose
    layout information is held in 'molinfo' */
ZMatrix::ZMatrix(const MoleculeInfoData &info)
        : ConcreteProperty<ZMatrix,MoleculeProperty>(),
          molinfo(info)
{}

/** Copy constructor */
ZMatrix::ZMatrix(const ZMatrix &other)
        : ConcreteProperty<ZMatrix,MoleculeProperty>(),
          molinfo(other.molinfo),
          zmat(other.zmat), atomidx_to_zmat(other.atomidx_to_zmat),
          zmat_build_order(other.zmat_build_order)
{}

/** Destructor */
ZMatrix::~ZMatrix()
{}

/** Copy assignment operator */
ZMatrix& ZMatrix::operator=(const ZMatrix &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        molinfo = other.molinfo;
        zmat = other.zmat;
        atomidx_to_zmat = other.atomidx_to_zmat;
        zmat_build_order = other.zmat_build_order;
    }

    return *this;
}

/** Comparison operator */
bool ZMatrix::operator==(const ZMatrix &other) const
{
    return molinfo == other.molinfo and zmat == other.zmat;
}

/** Comparison operator */
bool ZMatrix::operator!=(const ZMatrix &other) const
{
    return molinfo != other.molinfo or zmat != other.zmat;
}

/** Recalculate the optimum order with which to build atoms using
    this z-matrix - this order will ensure that atoms are built after
    any atoms on which they depend. This will raise an exception
    if a circular reference is detected

    \throw SireMol::zmatrix_error
*/
void ZMatrix::rebuildOrder()
{
    const ZMatrixLine *lines_array = zmat.constData();
    const int nlines = zmat.count();

    if (nlines == 0)
    {
        zmat_build_order = QVector<int>();
        return;
    }

    QVector<int> new_order( nlines );
    new_order.squeeze();

    QSet<int> unconstructed_atoms;
    unconstructed_atoms.reserve(nlines);

    for (int i=0; i<nlines; ++i)
    {
        unconstructed_atoms.insert( lines_array[i].atom().value() );
    }

    int nassigned = 0;
    int old_nassigned = 0;

    while (nassigned < nlines)
    {
        for (int i=0; i<nlines; ++i)
        {
            const ZMatrixLine &line = lines_array[i];

            if ( unconstructed_atoms.contains(line.atom()) and
                 not (unconstructed_atoms.contains(line.bond().value()) or
                      unconstructed_atoms.contains(line.angle().value()) or
                      unconstructed_atoms.contains(line.dihedral().value())) )
            {
                //this atom does not depend on any unconstructed atoms!
                new_order[nassigned] = i;
                ++nassigned;
                unconstructed_atoms.remove(line.atom().value());
            }
        }

        if (nassigned == old_nassigned)
        {
            //no atoms were assigned on this loop - this means that there
            //are circular dependencies
            throw SireMove::zmatrix_error( QObject::tr(
                    "The z-matrix contains a circular dependency (involving "
                    "the atoms with indicies %1.\n%2")
                        .arg(Sire::toString(unconstructed_atoms))
                        .arg(this->toString()), CODELOC );
        }
    }

    zmat_build_order = new_order;
}

/** Return the layout of the molecule whose z-matrix is contained
    in this object */
const MoleculeInfoData& ZMatrix::info() const
{
    if (molinfo.constData() == 0)
        return MoleculeInfoData::null();
    else
        return *molinfo;
}

/** Return a string representation of this z-matrix */
QString ZMatrix::toString() const
{
    QStringList lines;

    lines.append( QObject::tr("ZMatrix nAtoms() == %1").arg(info().nAtoms()) );

    for (int i=0; i<zmat.count(); ++i)
    {
        const ZMatrixLine &line = zmat.at(i);

        lines.append( QObject::tr("%1(%5)-%2(%6)-%3(%7)-%4(%8) : %9 %10 %11")
                        .arg( info().name(line.atom()).value(),
                              info().name(line.bond()).value(),
                              info().name(line.angle()).value(),
                              info().name(line.dihedral()).value() )
                        .arg(line.atom().value())
                        .arg(line.bond().value())
                        .arg(line.angle().value())
                        .arg(line.dihedral().value())
                        .arg(line.bondDelta().to(angstroms))
                        .arg(line.angleDelta().to(degrees))
                        .arg(line.dihedralDelta().to(degrees)) );
    }

    return lines.join("\n");
}

/** Return the line for the atom identified by 'atom'. This
    raises an exception if this is not a valid ID, or if
    there is no z-matrix line for this atom

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
const ZMatrixLine& ZMatrix::operator[](const AtomID &atom) const
{
    AtomIdx idx = info().atomIdx(atom);

    if (not atomidx_to_zmat.contains(idx))
        throw SireMove::zmatrix_error( QObject::tr(
            "The atom with ID %1 does not appear in the z-matrix.")
                .arg(atom.toString()), CODELOC );

    return zmat.at( atomidx_to_zmat.value(idx) );
}

/** Return the line for the atom identified by 'atom'. This
    raises an exception if this is not a valid ID, or if
    there is no z-matrix line for this atom

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
const ZMatrixLine& ZMatrix::at(const AtomID &atom) const
{
    return this->operator[](atom);
}

/** Return all of the lines of the z-matrix */
const QVector<ZMatrixLine>& ZMatrix::lines() const
{
    return zmat;
}

/** Return the index of AtomIdx to z-matrix line number. This
    is used to index the output of ZMatrix::lines() */
const QHash<AtomIdx,int>& ZMatrix::index() const
{
    return atomidx_to_zmat;
}

/** Return whether or not this z-matrix contains a line for the
    atom with ID 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const AtomID &atom) const
{
    return atomidx_to_zmat.contains( info().atomIdx(atom) );
}

static bool zmatrixDependsOnAtom(const QVector<ZMatrixLine> &zmat,
                                 const AtomIdx &atom)
{
    const ZMatrixLine *lines_array = zmat.constData();
    const int nlines = zmat.count();

    for (int i=0; i<nlines; ++i)
    {
        const ZMatrixLine &line = lines_array[i];

        if (atom == line.bond() or atom == line.angle() or
            atom == line.dihedral())
        {
            return true;
        }
    }

    return false;
}

/** Add a line to the z-matrix that gives the coordinates of the
    atom 'atom' based on the passed bond, angle and dihedal atoms

    An exception is raised if adding this line would lead
    to a circular reference

    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::add(const AtomID &atom, const AtomID &bond,
                  const AtomID &angle, const AtomID &dihedral)
{
    AtomIdx atm = info().atomIdx(atom);
    AtomIdx bnd = info().atomIdx(bond);
    AtomIdx ang = info().atomIdx(angle);
    AtomIdx dih = info().atomIdx(dihedral);

    if (atm == bnd or atm == ang or atm == dih or
        bnd == ang or bnd == dih or
        ang == dih)
    {
        throw SireMove::zmatrix_error( QObject::tr(
            "You cannot add a z-matrix line with repeated atoms! "
            "(%1-%2-%3-%4)")
                .arg(atom.toString(), bond.toString(),
                     angle.toString(), dihedral.toString()), CODELOC );
    }

    const int nlines = zmat.count();

    if (zmatrixDependsOnAtom(zmat, atm))
    {
        ZMatrix old_state(*this);

        try
        {
            zmat.append( ZMatrixLine(atm,bnd,ang,dih) );
            atomidx_to_zmat.insert(atm, nlines);
            this->rebuildOrder();
        }
        catch(...)
        {
            ZMatrix::operator=(old_state);
            throw;
        }
    }
    else
    {
        //we can just append this line onto the end of the z-matrix
        zmat.append( ZMatrixLine(atm,bnd,ang,dih) );
        atomidx_to_zmat.insert(atm, nlines);
        zmat_build_order.append(nlines);
    }
}

/** Add a line to the z-matrix that gives the coordinates of the
    passed dihedral

    An exception is raised if adding this line would lead
    to a circular reference

    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::add(const DihedralID &dihedral)
{
    this->add( dihedral.atom0(), dihedral.atom1(),
               dihedral.atom2(), dihedral.atom3() );
}

void ZMatrix::reindex()
{
    const int nlines = zmat.count();
    const ZMatrixLine *lines_array = zmat.constData();

    atomidx_to_zmat.clear();
    atomidx_to_zmat.reserve(nlines);

    for (int i=0; i<nlines; ++i)
    {
        atomidx_to_zmat.insert( lines_array[i].atom(), i );
    }

    this->rebuildOrder();
}

/** Remove the z-matrix line that gives the coordinates of the
    atom 'atom'. This does nothing if there is no z-matrix
    line for this atom - note this removes the lines for
    *all* matching atoms

    \throw SireMol::missing_atom
    \throw SireError::invalid_index
*/
void ZMatrix::remove(const AtomID &atom)
{
    QList<AtomIdx> atomidxs = atom.map(info());

    bool need_reindex = false;

    foreach (AtomIdx atomidx, atomidxs)
    {
        if (atomidx_to_zmat.contains(atomidx))
        {
            zmat.remove( atomidx_to_zmat.value(atomidx) );
            need_reindex = true;
        }
    }

    if (need_reindex)
        this->reindex();
}

/** Remove the z-matrix line involving the specified atom, bond,
    angle and dihedral. This matches only a single line, and
    will do nothing if this line is not in this z-matrix

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void ZMatrix::remove(const AtomID &atom, const AtomID &bond,
                     const AtomID &angle, const AtomID &dihedral)
{
    AtomIdx atm = info().atomIdx(atom);
    AtomIdx bnd = info().atomIdx(bond);
    AtomIdx ang = info().atomIdx(angle);
    AtomIdx dih = info().atomIdx(dihedral);

    if (atomidx_to_zmat.contains(atm))
    {
        const ZMatrixLine &line = zmat.at( atomidx_to_zmat.value(atm) );

        if (line.bond() == bnd and line.angle() == ang and
            line.dihedral() == dih)
        {
            zmat.remove( atomidx_to_zmat.value(atm) );
            this->reindex();
        }
    }
}

/** Remove the z-matrix line involving the specified dihedral.
    This matches only a single line, and
    will do nothing if this line is not in this z-matrix

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
void ZMatrix::remove(const DihedralID &dihedral)
{
    this->remove( dihedral.atom0(), dihedral.atom1(),
                  dihedral.atom2(), dihedral.atom3() );
}

/** Add the z-matrix line 'zmatline'.

    An exception is raised if adding this line would lead to
    a circular reference

    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::add(const ZMatrixLine &zmatline)
{
    this->add( zmatline.atom(), zmatline.bond(),
               zmatline.angle(), zmatline.dihedral() );

    int i = atomidx_to_zmat.value( info().atomIdx(zmatline.atom()) );

    zmat[i].setBondDelta( zmatline.bondDelta() );
    zmat[i].setAngleDelta( zmatline.angleDelta() );
    zmat[i].setDihedralDelta( zmatline.dihedralDelta() );
}

/** Remove the z-matrix line 'zmatline' from this z-matrix */
void ZMatrix::remove(const ZMatrixLine &zmatline)
{
    this->remove( zmatline.atom(), zmatline.bond(),
                  zmatline.angle(), zmatline.dihedral() );
}

/** Return the index of the z-matrix line that positions the atom
    with ID 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const AtomID &atom) const
{
    AtomIdx atomidx = info().atomIdx(atom);

    if (not atomidx_to_zmat.contains(atomidx))
    {
        throw SireMove::zmatrix_error( QObject::tr(
            "The atom %1 is not present in the z-matrix")
                .arg(atom.toString()),
                    CODELOC );
    }

    return atomidx_to_zmat.value(atomidx);
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const AtomID &atom, const AtomID &bond) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);

    if (atomidx_to_zmat.contains(atomidx))
    {
        if ( zmat.at(atomidx_to_zmat.value(atomidx)).bond() == bondidx )
            return atomidx_to_zmat.value(atomidx);
    }

    if (atomidx_to_zmat.contains(bondidx))
    {
        if ( zmat.at(atomidx_to_zmat.value(bondidx)).bond() == atomidx )
            return atomidx_to_zmat.value(bondidx);
    }

    throw SireMove::zmatrix_error( QObject::tr(
            "The bond %1-%2 is not present in the z-matrix.")
                .arg(atom.toString(), bond.toString()), CODELOC );

    return -1;
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'-'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const AtomID &atom, const AtomID &bond,
                      const AtomID &angle) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);
    AtomIdx angleidx = info().atomIdx(angle);

    if (atomidx_to_zmat.contains(atomidx))
    {
        int idx = atomidx_to_zmat.value(atomidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == angleidx)
            return idx;
    }

    if (atomidx_to_zmat.contains(angleidx))
    {
        int idx = atomidx_to_zmat.value(angleidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == atomidx)
            return idx;
    }

    throw SireMove::zmatrix_error( QObject::tr(
            "The angle %1-%2-%3 is not present in the z-matrix.")
                .arg(atom.toString(), bond.toString(),
                     angle.toString()), CODELOC );

    return -1;
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'-'angle'-'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const AtomID &atom, const AtomID &bond,
                      const AtomID &angle, const AtomID &dihedral) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);
    AtomIdx angleidx = info().atomIdx(angle);
    AtomIdx dihedralidx = info().atomIdx(dihedral);

    if (atomidx_to_zmat.contains(atomidx))
    {
        int idx = atomidx_to_zmat.value(atomidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == angleidx and
            line.dihedral() == dihedralidx)
            return idx;
    }

    if (atomidx_to_zmat.contains(dihedralidx))
    {
        int idx = atomidx_to_zmat.value(dihedralidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == angleidx and line.angle() == bondidx and
            line.dihedral() == atomidx)
            return idx;
    }

    throw SireMove::zmatrix_error( QObject::tr(
            "The dihedral %1-%2-%3-%4 is not present in the z-matrix.")
                .arg(atom.toString(), bond.toString(),
                     angle.toString(), dihedral.toString()), CODELOC );

    return -1;
}

/** Return whether or not this z-matrix contains the bond between
    atoms 'atom'-'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const AtomID &atom, const AtomID &bond) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);

    if (atomidx_to_zmat.contains(atomidx))
    {
        if ( zmat.at(atomidx_to_zmat.value(atomidx)).bond() == bondidx )
            return true;;
    }

    if (atomidx_to_zmat.contains(bondidx))
    {
        if ( zmat.at(atomidx_to_zmat.value(bondidx)).bond() == atomidx )
            return true;
    }

    return false;
}

/** Return whether or not this z-matrix contains the angle
    'atom'-'bond'-'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const AtomID &atom, const AtomID &bond, const AtomID &angle) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);
    AtomIdx angleidx = info().atomIdx(angle);

    if (atomidx_to_zmat.contains(atomidx))
    {
        int idx = atomidx_to_zmat.value(atomidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == angleidx)
            return true;
    }

    if (atomidx_to_zmat.contains(angleidx))
    {
        int idx = atomidx_to_zmat.value(angleidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == atomidx)
            return true;
    }

    return false;
}

/** Return whether or not this z-matrix contains the dihedral
    'atom'-'bond'-'angle'-'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const AtomID &atom, const AtomID &bond,
                       const AtomID &angle, const AtomID &dihedral) const
{
    AtomIdx atomidx = info().atomIdx(atom);
    AtomIdx bondidx = info().atomIdx(bond);
    AtomIdx angleidx = info().atomIdx(angle);
    AtomIdx dihedralidx = info().atomIdx(dihedral);

    if (atomidx_to_zmat.contains(atomidx))
    {
        int idx = atomidx_to_zmat.value(atomidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == bondidx and line.angle() == angleidx and
            line.dihedral() == dihedralidx)
            return true;
    }

    if (atomidx_to_zmat.contains(dihedralidx))
    {
        int idx = atomidx_to_zmat.value(dihedralidx);
        const ZMatrixLine &line = zmat.at(idx);

        if (line.bond() == angleidx and line.angle() == bondidx and
            line.dihedral() == atomidx)
            return true;
    }

    return false;
}

/** Return whether or not this z-matrix contains 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const BondID &bond) const
{
    return this->contains(bond.atom0(), bond.atom1());
}

/** Return whether or not this z-matrix contains 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const AngleID &angle) const
{
    return this->contains(angle.atom0(), angle.atom1(), angle.atom2());
}

/** Return whether or not this z-matrix contains 'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrix::contains(const DihedralID &dihedral) const
{
    return this->contains(dihedral.atom0(), dihedral.atom1(),
                          dihedral.atom2(), dihedral.atom3());
}

/** Return the index of the z-matrix line that defines the
    bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const BondID &bond) const
{
    return this->getIndex(bond.atom0(), bond.atom1());
}

/** Return the index of the z-matrix line that defines the
    angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const AngleID &angle) const
{
    return this->getIndex(angle.atom0(), angle.atom1(), angle.atom2());
}

/** Return the index of the z-matrix line that defines the
    dihedral 'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrix::getIndex(const DihedralID &dihedral) const
{
    return this->getIndex(dihedral.atom0(), dihedral.atom1(),
                          dihedral.atom2(), dihedral.atom3());
}

/** Set the maximum amount that the bond for the atom 'atom'
    can be moved to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setBondDelta(const AtomID &atom, const Length &delta)
{
    int idx = this->getIndex(atom);
    zmat[idx].setBondDelta(delta);
}

/** Set the maximum amount that the angle for the atom 'atom'
    can be changed to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setAngleDelta(const AtomID &atom, const Angle &delta)
{
    int idx = this->getIndex(atom);
    zmat[idx].setAngleDelta(delta);
}

/** Set the maximum amount that the dihedral for the atom 'atom'
    can be changed to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setDihedralDelta(const AtomID &atom, const Angle &delta)
{
    int idx = this->getIndex(atom);
    zmat[idx].setDihedralDelta(delta);
}

/** Set the maximum amount that the bond between atoms 'atom'-'bond'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setBondDelta(const AtomID &atom, const AtomID &bond,
                           const Length &delta)
{
    int idx = getIndex(atom, bond);
    zmat[idx].setBondDelta(delta);
}

/** Set the maximum amount that the angle between atoms 'atom'-'bond'-'angle'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setAngleDelta(const AtomID &atom, const AtomID &bond,
                            const AtomID &angle, const Angle &delta)
{
    int idx = getIndex(atom, bond, angle);
    zmat[idx].setAngleDelta(delta);
}

/** Set the maximum amount that the dihedral between atoms
    'atom'-'bond'-'angle'-'dihedral' can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setDihedralDelta(const AtomID &atom, const AtomID &bond,
                               const AtomID &angle, const AtomID &dihedral,
                               const Angle &delta)
{
    int idx = getIndex(atom, bond, angle, dihedral);
    zmat[idx].setDihedralDelta(delta);
}

/** Set the maximum amount that the bond 'bond'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setDelta(const BondID &bond, const Length &delta)
{
    this->setBondDelta( bond.atom0(), bond.atom1(), delta );
}

/** Set the maximum amount that the angle 'angle'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setDelta(const AngleID &angle, const Angle &delta)
{
    this->setAngleDelta( angle.atom0(), angle.atom1(), angle.atom2(), delta );
}

/** Set the maximum amount that the dihedral 'dihedral'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrix::setDelta(const DihedralID &dihedral, const Angle &delta)
{
    this->setDihedralDelta( dihedral.atom0(), dihedral.atom1(),
                            dihedral.atom2(), dihedral.atom3(), delta );
}

/** Return the maximum amount that the bond to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrix::bondDelta(const AtomID &atom) const
{
    return zmat[getIndex(atom)].bondDelta();
}

/** Return the maximum amount that the angle to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::angleDelta(const AtomID &atom) const
{
    return zmat[getIndex(atom)].angleDelta();
}

/** Return the maximum amount that the dihedral to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::dihedralDelta(const AtomID &atom) const
{
    return zmat[getIndex(atom)].dihedralDelta();
}

/** Return the maximum amount that the bond between atoms
    'atom'-'bond' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrix::bondDelta(const AtomID &atom, const AtomID &bond) const
{
    return zmat[getIndex(atom,bond)].bondDelta();
}

/** Return the maximum amount that the angle between atoms
    'atom'-'bond'-'angle' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::angleDelta(const AtomID &atom, const AtomID &bond,
                          const AtomID &angle) const
{
    return zmat[getIndex(atom,bond,angle)].angleDelta();
}

/** Return the maximum amount that the dihedral between atoms
    'atom'-'bond'-'angle'-'dihedral' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::dihedralDelta(const AtomID &atom, const AtomID &bond,
                             const AtomID &angle, const AtomID &dihedral) const
{
    return zmat[getIndex(atom,bond,angle,dihedral)].dihedralDelta();
}

/** Return the maximum amount that 'bond' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrix::delta(const BondID &bond) const
{
    return this->bondDelta( bond.atom0(), bond.atom1() );
}

/** Return the maximum amount that 'angle' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::delta(const AngleID &angle) const
{
    return this->angleDelta( angle.atom0(), angle.atom1(), angle.atom2() );
}

/** Return the maximum amount that 'dihedral' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrix::delta(const DihedralID &dihedral) const
{
    return this->dihedralDelta( dihedral.atom0(), dihedral.atom1(),
                                dihedral.atom2(), dihedral.atom3() );
}

/** Return whether or not this z-matrix is compatible with the
    the molecule whose info is in 'molinfo' */
bool ZMatrix::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return info().UID() == molinfo.UID();
}

/** Return a z-matrix that only contains lines that involve the atoms
    that are in 'selection' */
ZMatrix ZMatrix::matchToSelection(const AtomSelection &selection) const
{
    selection.assertCompatibleWith(*molinfo);

    if (selection.selectedAll())
    {
        ZMatrix zmat(*this);
        zmat.molinfo = selection.info();
        return zmat;
    }
    else if (selection.selectedNone())
    {
        ZMatrix zmatrix;
        zmatrix.molinfo = selection.info();
        return zmatrix;
    }

    QVector<AtomIdx> selected_atoms = selection.selectedAtoms();

    int nselected = selected_atoms.count();

    QVector<int> selected_atom_lines;
    selected_atom_lines.reserve(nselected);

    const AtomIdx *selected_atoms_array = selected_atoms.constData();

    for (int i=0; i<nselected; ++i)
    {
        if (atomidx_to_zmat.contains(selected_atoms_array[i]))
        {
            selected_atom_lines.append(
                              atomidx_to_zmat.value(selected_atoms_array[i]) );
        }
    }

    qSort(selected_atom_lines.begin(), selected_atom_lines.end());

    int nzmat = selected_atom_lines.count();
    const int *selected_atom_lines_array = selected_atom_lines.constData();

    QVector<ZMatrixLine> new_zmat(nzmat);
    ZMatrixLine *new_zmat_array = new_zmat.data();
    const ZMatrixLine *old_zmat_array = zmat.constData();

    QHash<AtomIdx,int> new_index;
    new_index.reserve(nzmat);

    for (int i=0; i<nzmat; ++i)
    {
        new_zmat_array[i] = old_zmat_array[ selected_atom_lines_array[i] ];
        new_index.insert( new_zmat_array[i].atom(), i );
    }

    ZMatrix zmatrix(*this);
    zmatrix.zmat = new_zmat;
    zmatrix.molinfo = selection.info();
    zmatrix.atomidx_to_zmat = new_index;
    zmatrix.rebuildOrder();

    return zmatrix;
}

/** Return the number of lines in this z-matrix */
int ZMatrix::nLines() const
{
    return zmat.count();
}

/** Return the correct atom build order for this z-matrix */
const QVector<int>& ZMatrix::atomBuildOrder() const
{
    return zmat_build_order;
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom matching
    functions in 'atommatcher' to match atoms from the current molecule
    to the atoms in the molecule whose layout is in 'molinfo'

    This only copies the z-matrix lines when all of the atoms in that
    line have been matched. It does not copy lines where there is no
    match. Use the ZMatrix::nLines() function to check that the number
    of lines doesn't change if you want to ensure that all of the
    lines have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr ZMatrix::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                             const AtomMatcher &atommatcher) const
{
    if (not atommatcher.changesOrder(this->info(), molinfo))
    {
        //the order and number of atoms is the same - the z-matrix
        //can just be copied
        ZMatrix ret(molinfo);

        ret.zmat = zmat;
        ret.atomidx_to_zmat = atomidx_to_zmat;

        return ret;
    }

    QHash<AtomIdx,AtomIdx> matched_atoms = atommatcher.match(this->info(), molinfo);

    return this->_pvt_makeCompatibleWith(molinfo, matched_atoms);
}

/** Return a copy of this property that has been made to be compatible
    with the molecule layout in 'molinfo' - this uses the atom mapping
    in 'map' to match atoms from the current molecule to the atoms in
    the molecule whose layout is in 'molinfo'

    This only copies the z-matrix lines when all of the atoms in that
    line have been matched. It does not copy lines where there is no
    match. Use the ZMatrix::nLines() function to check that the number
    of lines doesn't change if you want to ensure that all of the
    lines have been copied.

    \throw SireError::incompatible_error
*/
PropertyPtr ZMatrix::_pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                             const QHash<AtomIdx,AtomIdx> &map) const
{
    ZMatrix ret(molinfo);

    int nlines = zmat.count();

    for (int i=0; i<nlines; ++i)
    {
        const ZMatrixLine &line = zmat.at(i);

        AtomIdx atm = map.value(line.atom(), AtomIdx(-1));
        AtomIdx bnd = map.value(line.bond(), AtomIdx(-1));
        AtomIdx ang = map.value(line.angle(), AtomIdx(-1));
        AtomIdx dih = map.value(line.dihedral(), AtomIdx(-1));

        if (atm == -1 or bnd == -1 or ang == -1 or dih == -1)
            continue;

        ret.add( atm, bnd, ang, dih );
        ret.setBondDelta( atm, line.bondDelta() );
        ret.setAngleDelta( atm, line.angleDelta() );
        ret.setDihedralDelta( atm, line.dihedralDelta() );
    }

    return ret;
}

const char* ZMatrix::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ZMatrix>() );
}

//////////
////////// Implementation of ZMatrixCoords
//////////

static const RegisterMetaType<ZMatrixCoords> r_zmatcoords;

/** Serialise to a binary datastream */
QDataStream SIREMOVE_EXPORT &operator<<(QDataStream &ds, const ZMatrixCoords &zmatcoords)
{
    writeHeader(ds, r_zmatcoords, 1);

    SharedDataStream sds(ds);

    sds << zmatcoords.zmat << zmatcoords.internal_coords
        << zmatcoords.cartesian_coords
        << static_cast<const MoleculeProperty&>(zmatcoords);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMOVE_EXPORT &operator>>(QDataStream &ds, ZMatrixCoords &zmatcoords)
{
    VersionID v = readHeader(ds, r_zmatline);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> zmatcoords.zmat >> zmatcoords.internal_coords
            >> zmatcoords.cartesian_coords
            >> static_cast<MoleculeProperty&>(zmatcoords);
    }
    else
        throw version_error( v, "1", r_zmatcoords, CODELOC );

    return ds;
}

/** Null constructor */
ZMatrixCoords::ZMatrixCoords()
              : ConcreteProperty<ZMatrixCoords,MoleculeProperty>(),
                need_rebuild(false)
{}

/** Internal function used to calculate the internal coordinates
    of the passed z-matrix line */
Vector ZMatrixCoords::getInternalCoords(const ZMatrixLine &line) const
{
    Vector atom = cartesian_coords[info().cgAtomIdx(line.atom())];
    Vector bond = cartesian_coords[info().cgAtomIdx(line.bond())];
    Vector angle = cartesian_coords[info().cgAtomIdx(line.angle())];
    Vector dihedral = cartesian_coords[info().cgAtomIdx(line.dihedral())];

    return Vector( Vector::distance(atom,bond),
                   Vector::angle(atom,bond,angle).value(),
                   Vector::dihedral(atom,bond,angle,dihedral).value() );
}

/** Internal function used to rebuild the internal coordinates
    from the cartesian coordinates */
void ZMatrixCoords::rebuildInternals()
{
    if (zmat.lines().isEmpty())
    {
        internal_coords.clear();
        return;
    }

    this->rebuildCartesian();

    int nlines = zmat.lines().count();

    const ZMatrixLine *lines_array = zmat.lines().constData();

    internal_coords = QVector<Vector>(nlines);
    internal_coords.squeeze();

    Vector *internal_coords_array = internal_coords.data();

    for (int i=0; i<nlines; ++i)
    {
        internal_coords_array[i] = this->getInternalCoords(lines_array[i]);
    }
}

/** Construct the z-matrix for the molecule 'molecule' using
    the passed property map to find the coordinates property, and
    also picking up the molecule's existing z-matrix */
ZMatrixCoords::ZMatrixCoords(const PartialMolecule &molecule, const PropertyMap &map)
              : ConcreteProperty<ZMatrixCoords,MoleculeProperty>(),
                zmat(molecule.molecule()),
                need_rebuild(false)
{
    QTime t;
    t.start();

    cartesian_coords = molecule.molecule().property( map["coordinates"] )
                                          .asA<AtomCoords>();

    const PropertyName &zmatrix_property = map["z-matrix"];

    if (molecule.hasProperty(zmatrix_property))
    {
        zmat = molecule.molecule().property(zmatrix_property).asA<ZMatrix>();

        if (not molecule.selection().selectedAll())
            zmat = zmat.matchToSelection( molecule.selection() );

        //calculate the internal coordinates
        this->rebuildInternals();
    }
}

/** Construct the z-matrix for the passed coordinates and z-matrix */
ZMatrixCoords::ZMatrixCoords(const ZMatrix &zmatrix, const AtomCoords &coords)
              : ConcreteProperty<ZMatrixCoords,MoleculeProperty>(),
                zmat(zmatrix), cartesian_coords(coords), need_rebuild(false)
{
    coords.assertCompatibleWith( zmat.info() );

    //calculate the internal coordinates
    this->rebuildInternals();
}

/** Construct the z-matrix for the molecule 'molecule' using the
    z-matrix in 'zmatrix' and using the passed property map to
    find the coordinates */
ZMatrixCoords::ZMatrixCoords(const ZMatrix &zmatrix, const PartialMolecule &molecule,
                             const PropertyMap &map)
              : ConcreteProperty<ZMatrixCoords,MoleculeProperty>(),
                zmat(zmatrix), need_rebuild(false)
{
    zmatrix.assertCompatibleWith( molecule.data().info() );

    cartesian_coords = molecule.molecule().property( map["coordinates"] )
                                          .asA<AtomCoords>();

    if (not molecule.selection().selectedAll())
        zmat = zmatrix.matchToSelection( molecule.selection() );

    //calculate the internal coordinates
    this->rebuildInternals();
}

/** Copy constructor */
ZMatrixCoords::ZMatrixCoords(const ZMatrixCoords &other)
              : ConcreteProperty<ZMatrixCoords,MoleculeProperty>(other),
                zmat(other.zmat), internal_coords(other.internal_coords),
                cartesian_coords(other.cartesian_coords),
                need_rebuild(other.need_rebuild)
{}

/** Destructor */
ZMatrixCoords::~ZMatrixCoords()
{}

/** Copy assignment operator */
ZMatrixCoords& ZMatrixCoords::operator=(const ZMatrixCoords &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        zmat = other.zmat;
        internal_coords = other.internal_coords;
        cartesian_coords = other.cartesian_coords;
        need_rebuild = other.need_rebuild;
    }

    return *this;
}

/** Comparison operator */
bool ZMatrixCoords::operator==(const ZMatrixCoords &other) const
{
    return this == &other or
           (zmat == other.zmat and internal_coords == other.internal_coords and
            cartesian_coords == other.cartesian_coords and
            need_rebuild == other.need_rebuild);
}

/** Comparison operator */
bool ZMatrixCoords::operator!=(const ZMatrixCoords &other) const
{
    return not this->operator==(other);
}

Q_GLOBAL_STATIC( QMutex, zmatrixMutex );

/** Internal function called to rebuild the cartesian coordinates */
void ZMatrixCoords::_pvt_rebuildCartesian()
{
    int nlines = zmat.lines().count();

    const ZMatrixLine *lines_array = zmat.lines().constData();
    const int *build_order = zmat.atomBuildOrder().constData();
    const Vector *internal_coords_array = internal_coords.constData();

    BOOST_ASSERT( zmat.atomBuildOrder().count() == nlines );

    AtomCoords new_coords = cartesian_coords;

    for (int i=0; i<nlines; ++i)
    {
        int build_atom = build_order[i];

        const ZMatrixLine &line = lines_array[build_atom];
        const Vector &internal = internal_coords_array[build_atom];

        //get the coordinates of the bond, angle and dihedral atoms
        Vector bond = new_coords[ info().cgAtomIdx(line.bond()) ];
        Vector angle = new_coords[ info().cgAtomIdx(line.angle()) ];
        Vector dihedral = new_coords[ info().cgAtomIdx(line.dihedral()) ];

        //now use these to build the coordinates of the atom
        new_coords.set( info().cgAtomIdx(line.atom()),
                        Vector::generate(internal[0], bond,
                                         Angle(internal[1]), angle,
                                         Angle(internal[2]), dihedral) );
    }

    cartesian_coords = new_coords;
}

/** Internal function called to rebuild the cartesian coordinates
    from the internal coordinates */
void ZMatrixCoords::rebuildCartesian() const
{
    //need a mutex as we have declared this as a const function
    //and this breaks implicit sharing
    QMutexLocker lkr( zmatrixMutex() );

    if (not need_rebuild)
        return;

    const_cast<ZMatrixCoords*>(this)->_pvt_rebuildCartesian();
}

/** Return the z-matrix line for the atom identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
ZMatrixCoordsLine ZMatrixCoords::operator[](const AtomID &atom) const
{
    AtomIdx idx = zmat.info().atomIdx(atom);

    if (not zmat.index().contains(idx))
        throw SireMove::zmatrix_error( QObject::tr(
            "The atom %1 does not appear in the z-matrix.")
                .arg(atom.toString()), CODELOC );

    const Vector &coords = internal_coords.at( zmat.index().value(idx) );

    return ZMatrixCoordsLine( zmat[idx], Length(coords[0]),
                              Angle(coords[1]), Angle(coords[2]) );
}

/** Return the z-matrix line for the atom identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
ZMatrixCoordsLine ZMatrixCoords::at(const AtomID &atom) const
{
    return this->operator[](atom);
}

/** Return all of the lines of the z-matrix */
QVector<ZMatrixCoordsLine> ZMatrixCoords::lines() const
{
    if (zmat.lines().isEmpty())
        return QVector<ZMatrixCoordsLine>();

    QVector<ZMatrixCoordsLine> coords( zmat.lines().count() );

    const ZMatrixLine *zmat_array = zmat.lines().constData();
    const Vector *internal_coords_array = internal_coords.constData();
    ZMatrixCoordsLine *coords_array = coords.data();

    for (int i=0; i<zmat.lines().count(); ++i)
    {
        const Vector &c = internal_coords_array[i];

        coords_array[i] = ZMatrixCoordsLine( zmat_array[i], Length(c[0]),
                                             Angle(c[1]), Angle(c[2]) );
    }

    return coords;
}

/** Return just the internal coordinates */
const QVector<Vector>& ZMatrixCoords::internalCoordinates() const
{
    return internal_coords;
}

/** Return the index of AtomIdx to z-matrix line number. This
    is used to index the output of ZMatrixCoords::lines() */
const QHash<AtomIdx,int>& ZMatrixCoords::index() const
{
    return zmat.index();
}

/** Return the index of the z-matrix line that positions the atom
    with ID 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const AtomID &atom) const
{
    return zmat.getIndex(atom);
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const AtomID &atom, const AtomID &bond) const
{
    return zmat.getIndex(atom, bond);
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'-'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const AtomID &atom, const AtomID &bond,
                            const AtomID &angle) const
{
    return zmat.getIndex(atom, bond, angle);
}

/** Return the index of the z-matrix line that positions the
    bond between atoms 'atom'-'bond'-'angle'-'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const AtomID &atom, const AtomID &bond,
                            const AtomID &angle, const AtomID &dihedral) const
{
    return zmat.getIndex(atom, bond, angle, dihedral);
}

/** Return the index of the z-matrix line that
    defines the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const BondID &bond) const
{
    return this->getIndex( bond.atom0(), bond.atom1() );
}

/** Return the index of the z-matrix line that
    defines the angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const AngleID &angle) const
{
    return this->getIndex( angle.atom0(), angle.atom1(), angle.atom2() );
}

/** Return the index of the z-matrix line that
    defines the dihedral 'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
int ZMatrixCoords::getIndex(const DihedralID &dihedral) const
{
    return this->getIndex( dihedral.atom0(), dihedral.atom1(),
                           dihedral.atom2(), dihedral.atom3() );
}

/** Return the layout of the molecule whose z-matrix is contained
    in this object */
const MoleculeInfoData& ZMatrixCoords::info() const
{
    return zmat.info();
}

/** Return the raw z-matrix, which does not contain any
    coordinate information */
const ZMatrix& ZMatrixCoords::zmatrix() const
{
    return zmat;
}

/** Return a string representation of this z-matrix */
QString ZMatrixCoords::toString() const
{
    QStringList output;

    int nats = info().nAtoms();

    output.append( QObject::tr("ZMatrix nAtoms() == %1").arg(nats) );

    QVector<ZMatrixCoordsLine> zmatrix = this->lines();

    for (AtomIdx i(0); i<nats; ++i)
    {
        if (zmat.index().contains(i))
        {
            const ZMatrixCoordsLine &line = zmatrix.at(zmat.index().value(i));

            output.append( QObject::tr("%1-%2-%3-%4 - %5 A : %6' : %7'")
                            .arg(line.atom().value())
                            .arg(line.bond().value())
                            .arg(line.angle().value())
                            .arg(line.dihedral().value())
                            .arg(line.bondLength())
                            .arg( line.angleSize().to(degrees) )
                            .arg( line.dihedralSize().to(degrees) ) );
        }
        else
        {
            output.append( QObject::tr("%1 - %2")
                                .arg(i)
                                .arg(cartesian_coords[info().cgAtomIdx(i)].toString()) );
        }
    }

    return output.join("\n");
}

/** Return whether or not the z-matrix contains an atom
    identified by 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const AtomID &atom) const
{
    return zmat.contains(atom);
}

/** Return whether or not the z-matrix defines
    the bond 'atom'-'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const AtomID &atom, const AtomID &bond) const
{
    return zmat.contains(atom, bond);
}

/** Return whether or not the z-matrix defines
    the angle 'atom'-'bond'-'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const AtomID &atom, const AtomID &bond,
                             const AtomID &angle) const
{
    return zmat.contains(atom, bond, angle);
}

/** Return whether or not the z-matrix defines
    the dihedral 'atom'-'bond'-'angle'-'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const AtomID &atom, const AtomID &bond,
                             const AtomID &angle, const AtomID &dihedral) const
{
    return zmat.contains(atom, bond, angle, dihedral);
}

/** Return whether or not the z-matrix defines
    the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const BondID &bond) const
{
    return zmat.contains(bond);
}

/** Return whether or not the z-matrix defines
    the angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const AngleID &angle) const
{
    return zmat.contains(angle);
}

/** Return whether or not the z-matrix defines
    the dihedral 'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
*/
bool ZMatrixCoords::contains(const DihedralID &dihedral) const
{
    return zmat.contains(dihedral);
}

/** Internal function used to add a single set of internal
    coordinates to the z-matrix */
void ZMatrixCoords::addInternal(const AtomIdx &atom)
{
    if (internal_coords.count() != zmat.index().count())
    {
        //this has added the coordinates
        if (internal_coords.count() != zmat.index().count() - 1)
            //more than one line has been added???
            throw SireError::program_bug( QObject::tr(
                "How can we have %1 lines of zmatrix and %2 lines of coordinates?")
                    .arg(zmat.index().count()).arg(internal_coords.count()),
                        CODELOC );

        //add space for the coordinates
        if (not zmat.index().contains(atom))
            throw SireError::program_bug( QObject::tr(
                "What's happened now? The index doesn't contain the atom %1.")
                    .arg(atom.toString()), CODELOC );

        internal_coords.insert( zmat.index().value(atom), 1, Vector() );
    }

    if (not zmat.index().contains(atom))
        throw SireError::program_bug( QObject::tr(
             "What's happened here? The index doesn't contain the atom %1.")
                 .arg(atom.toString()), CODELOC );

    //now calculate the coordinates
    int idx = zmat.index().value(atom);

    const ZMatrixLine &line = zmat.lines().at(idx);

    this->rebuildCartesian();

    internal_coords[idx] = this->getInternalCoords(line);
}

/** Add the z-matrix line containing the passed atoms. The internal
    coordinates of this line are set from the current cartesian coordinates

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index

    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::add(const AtomID &atom, const AtomID &bond,
                        const AtomID &angle, const AtomID &dihedral)
{
    ZMatrix old_zmat = zmat;

    try
    {
        zmat.add(atom, bond, angle, dihedral);
        this->addInternal( info().atomIdx(atom) );
    }
    catch(...)
    {
        zmat = old_zmat;
        throw;
    }
}

/** Add the z-matrix line containing the passed dihedral. The internal
    coordinates of this line are set from the current cartesian coordinates

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index

    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::add(const DihedralID &dihedral)
{
    this->add( dihedral.atom0(), dihedral.atom1(),
               dihedral.atom2(), dihedral.atom3() );
}

/** Add the z-matrix line using the supplied atoms and internal coordinates

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::add(const AtomID &atom,
                        const Length &bondlength, const AtomID &bond,
                        const Angle &anglesize, const AtomID &angle,
                        const Angle &dihedralsize, const AtomID &dihedral)
{
    ZMatrix old_zmat = zmat;

    try
    {
        zmat.add(atom, bond, angle, dihedral);

        AtomIdx atomidx = info().atomIdx(atom);

        this->addInternal(atomidx);

        BOOST_ASSERT( zmat.index().contains(atomidx) );

        int idx = zmat.index().value(atomidx);

        internal_coords[idx] = Vector( bondlength.value(), anglesize.value(),
                                       dihedralsize.value() );

        need_rebuild = true;
    }
    catch(...)
    {
        zmat = old_zmat;
        throw;
    }
}

/** Remove the z-matrix lines that build the atom(s) with ID 'atom' */
void ZMatrixCoords::remove(const AtomID &atom)
{
    this->rebuildCartesian();
    zmat.remove(atom);

    if (internal_coords.count() != zmat.index().count())
    {
        //something has been removed
        this->rebuildInternals();
    }
}

/** Remove the z-matrix line that builds 'atom' from 'bond',
    'angle' and 'dihedral' */
void ZMatrixCoords::remove(const AtomID &atom, const AtomID &bond,
                           const AtomID &angle, const AtomID &dihedral)
{
    this->rebuildCartesian();
    zmat.remove(atom, bond, angle, dihedral);

    if (internal_coords.count() != zmat.index().count())
    {
        this->rebuildInternals();
    }
}

/** Remove the z-matrix line that builds 'dihedral' */
void ZMatrixCoords::remove(const DihedralID &dihedral)
{
    this->remove( dihedral.atom0(), dihedral.atom1(),
                  dihedral.atom2(), dihedral.atom3() );
}

/** Add the z-matrix line 'zmatline'

    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::add(const ZMatrixLine &zmatline)
{
    ZMatrix old_zmat = zmat;

    try
    {
        zmat.add(zmatline);
        this->addInternal(zmatline.atom());
    }
    catch(...)
    {
        zmat = old_zmat;
        throw;
    }
}

/** Remove the z-matrix line 'zmatline' from this z-matrix */
void ZMatrixCoords::remove(const ZMatrixLine &zmatline)
{
    this->rebuildCartesian();
    zmat.remove(zmatline);

    if (internal_coords.count() != zmat.index().count())
    {
        //a line has been removed
        this->rebuildInternals();
    }
}

/** Add the z-matrix line 'zmatline' to this z-matrix, adding both
    the atoms and also setting the values of the internal coordinates

    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::add(const ZMatrixCoordsLine &zmatline)
{
    ZMatrix old_zmat = zmat;

    try
    {
        zmat.add(zmatline);

        this->addInternal(zmatline.atom());

        BOOST_ASSERT( zmat.index().contains(zmatline.atom()) );

        int idx = zmat.index().value(zmatline.atom());

        internal_coords[idx] = Vector( zmatline.bondLength().value(),
                                       zmatline.angleSize().value(),
                                       zmatline.dihedralSize().value() );

        need_rebuild = true;
    }
    catch(...)
    {
        zmat = old_zmat;
        throw;
    }
}

/** Return whether or not this zmatrix is compatible with the molecule
    with info 'molinfo' */
bool ZMatrixCoords::isCompatibleWith(const SireMol::MoleculeInfoData &molinfo) const
{
    return zmat.isCompatibleWith(molinfo);
}

/** Return the cartesian representation of these internal co-ordinates
    (convert from z-matrix coordinates to cartesian coordinates) */
const AtomCoords& ZMatrixCoords::toCartesian() const
{
    this->rebuildCartesian();
    return cartesian_coords;
}

/** Move the bond to the atom 'atom' by the length 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveBond(const AtomID &atom, const Length &delta)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setX( internal_coords[idx].x() + delta.value() );
    need_rebuild = true;
}

/** Move the angle to the atom 'atom' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveAngle(const AtomID &atom, const Angle &delta)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setY( internal_coords[idx].y() + delta.value() );
    need_rebuild = true;
}

/** Move the dihedral to the atom 'atom' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveDihedral(const AtomID &atom, const Angle &delta)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setZ( internal_coords[idx].z() + delta.value() );
    need_rebuild = true;
}

/** Change the bond between atoms 'atom0'-'atom1' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveBond(const AtomID &atom0, const AtomID &atom1,
                             const Length &delta)
{
    int idx = zmat.getIndex(atom0, atom1);
    internal_coords[idx].setX( internal_coords[idx].x() + delta.value() );
    need_rebuild = true;
}

/** Change the angle between atoms 'atom0'-'atom1'-'atom2' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveAngle(const AtomID &atom0, const AtomID &atom1,
                              const AtomID &atom2, const Angle &delta)
{
    int idx = zmat.getIndex(atom0, atom1, atom2);
    internal_coords[idx].setY( internal_coords[idx].y() + delta.value() );
    need_rebuild = true;
}

/** Change the dihedral between atoms 'atom0'-'atom1'-'atom2'-'atom3' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::moveDihedral(const AtomID &atom0, const AtomID &atom1,
                                 const AtomID &atom2, const AtomID &atom3,
                                 const Angle &delta)
{
    int idx = zmat.getIndex(atom0, atom1, atom2, atom3);
    internal_coords[idx].setZ( internal_coords[idx].z() + delta.value() );
    need_rebuild = true;
}

/** Change the bond 'bond' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::move(const BondID &bond, const Length &delta)
{
    this->moveBond( bond.atom0(), bond.atom1(), delta );
}

/** Change the angle 'angle' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::move(const AngleID &angle, const Angle &delta)
{
    this->moveAngle( angle.atom0(), angle.atom1(),
                     angle.atom2(), delta );
}

/** Change the dihedral 'dihedral' by 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::move(const DihedralID &dihedral, const Angle &delta)
{
    this->moveDihedral( dihedral.atom0(), dihedral.atom1(),
                        dihedral.atom2(), dihedral.atom3(), delta );
}

/** Set the bond to atom 'atom' to 'length'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setBond(const AtomID &atom, const Length &length)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setX( length.value() );
    need_rebuild = true;
}

/** Set the angle to atom 'atom' to 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setAngle(const AtomID &atom, const Angle &size)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setY( size.value() );
    need_rebuild = true;
}

/** Set the dihedral to atom 'atom' to 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDihedral(const AtomID &atom, const Angle &size)
{
    int idx = zmat.getIndex(atom);
    internal_coords[idx].setZ( size.value() );
    need_rebuild = true;
}

/** Set the bond between atoms 'atom0'-'atom1' to have
    the length 'length'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setBond(const AtomID &atom0, const AtomID &atom1,
                            const Length &length)
{
    int idx = zmat.getIndex(atom0, atom1);
    internal_coords[idx].setX( length.value() );
    need_rebuild = true;
}

/** Set the angle between atoms 'atom0'-'atom1'-'atom2' to have
    the size 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setAngle(const AtomID &atom0, const AtomID &atom1,
                             const AtomID &atom2, const Angle &size)
{
    int idx = zmat.getIndex(atom0, atom1, atom2);
    internal_coords[idx].setY( size.value() );
    need_rebuild = true;
}

/** Set the dihedral between atoms 'atom0'-'atom1'-'atom2'-'atom3' to have
    the size 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDihedral(const AtomID &atom0, const AtomID &atom1,
                                const AtomID &atom2, const AtomID &atom3,
                                const Angle &size)
{
    int idx = zmat.getIndex(atom0, atom1, atom2, atom3);
    internal_coords[idx].setZ( size.value() );
    need_rebuild = true;
}

/** Set the bond 'bond' to have the length 'length'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::set(const BondID &bond, const Length &length)
{
    this->setBond( bond.atom0(), bond.atom1(), length );
}

/** Set the angle 'angle' to have the size 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::set(const AngleID &angle, const Angle &size)
{
    this->setAngle( angle.atom0(), angle.atom1(), angle.atom2(), size );
}

/** Set the dihedral 'dihedral' to have the size 'size'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::set(const DihedralID &dihedral, const Angle &size)
{
    this->setDihedral( dihedral.atom0(), dihedral.atom1(),
                       dihedral.atom2(), dihedral.atom3(), size );
}

/** Return the length of the bond to the atom 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::bondLength(const AtomID &atom) const
{
    return Length(internal_coords[zmat.getIndex(atom)].x());
}

/** Return the size of the angle to the atom 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::angleSize(const AtomID &atom) const
{
    return Angle(internal_coords[zmat.getIndex(atom)].y());
}

/** Return the size of the dihedral to the atom 'atom'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::dihedralSize(const AtomID &atom) const
{
    return Angle(internal_coords[zmat.getIndex(atom)].z());
}

/** Return the length of the bond between the atoms in the z-matrix
    'atom'-'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::bondLength(const AtomID &atom, const AtomID &bond) const
{
    return Length(internal_coords[zmat.getIndex(atom,bond)].x());
}

/** Return the size of the angle between the atoms in the z-matrix
    'atom'-'bond'-'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::angleSize(const AtomID &atom, const AtomID &bond,
                               const AtomID &angle) const
{
    return Angle(internal_coords[zmat.getIndex(atom,bond,angle)].y());
}

/** Return the size of the dihedral between the atoms in the z-matrix
    'atom'-'bond'-'angle'-'dihedral'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::dihedralSize(const AtomID &atom, const AtomID &bond,
                                  const AtomID &angle, const AtomID &dihedral) const
{
    return Angle(internal_coords[zmat.getIndex(atom,bond,angle,dihedral)].z());
}

/** Return the length of the bond 'bond'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::length(const BondID &bond) const
{
    return this->bondLength( bond.atom0(), bond.atom1() );
}

/** Return the size of the angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::size(const AngleID &angle) const
{
    return this->angleSize( angle.atom0(), angle.atom1(), angle.atom2() );
}

/** Return the size of the angle 'angle'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::size(const DihedralID &dihedral) const
{
    return this->dihedralSize( dihedral.atom0(), dihedral.atom1(),
                               dihedral.atom2(), dihedral.atom3() );
}

/** Set the maximum amount that the bond for the atom 'atom'
    can be moved to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setBondDelta(const AtomID &atom, const Length &delta)
{
    zmat.setBondDelta(atom, delta);
}

/** Set the maximum amount that the angle for the atom 'atom'
    can be changed to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setAngleDelta(const AtomID &atom, const Angle &delta)
{
    zmat.setAngleDelta(atom, delta);
}

/** Set the maximum amount that the dihedral for the atom 'atom'
    can be changed to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDihedralDelta(const AtomID &atom, const Angle &delta)
{
    zmat.setDihedralDelta(atom, delta);
}

/** Set the maximum amount that the bond between atoms 'atom'-'bond'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setBondDelta(const AtomID &atom, const AtomID &bond,
                                 const Length &delta)
{
    zmat.setBondDelta(atom, bond, delta);
}

/** Set the maximum amount that the angle between atoms 'atom'-'bond'-'angle'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setAngleDelta(const AtomID &atom, const AtomID &bond,
                                  const AtomID &angle, const Angle &delta)
{
    zmat.setAngleDelta(atom, bond, angle, delta);
}

/** Set the maximum amount that the dihedral between atoms
    'atom'-'bond'-'angle'-'dihedral' can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDihedralDelta(const AtomID &atom, const AtomID &bond,
                                     const AtomID &angle, const AtomID &dihedral,
                                     const Angle &delta)
{
    zmat.setDihedralDelta(atom, bond, angle, dihedral, delta);
}

/** Set the maximum amount that the bond 'bond'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDelta(const BondID &bond, const Length &delta)
{
    this->setBondDelta( bond.atom0(), bond.atom1(), delta );
}

/** Set the maximum amount that the angle 'angle'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDelta(const AngleID &angle, const Angle &delta)
{
    this->setAngleDelta( angle.atom0(), angle.atom1(), angle.atom2(), delta );
}

/** Set the maximum amount that the dihedral 'dihedral'
    can be changed by to 'delta'

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
void ZMatrixCoords::setDelta(const DihedralID &dihedral, const Angle &delta)
{
    this->setDihedralDelta( dihedral.atom0(), dihedral.atom1(),
                            dihedral.atom2(), dihedral.atom3(), delta );
}

/** Return the maximum amount that the bond to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::bondDelta(const AtomID &atom) const
{
    return zmat.bondDelta(atom);
}

/** Return the maximum amount that the angle to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::angleDelta(const AtomID &atom) const
{
    return zmat.angleDelta(atom);
}

/** Return the maximum amount that the dihedral to atom 'atom'
    should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::dihedralDelta(const AtomID &atom) const
{
    return zmat.dihedralDelta(atom);
}

/** Return the maximum amount that the bond between atoms
    'atom'-'bond' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::bondDelta(const AtomID &atom, const AtomID &bond) const
{
    return zmat.bondDelta(atom, bond);
}

/** Return the maximum amount that the angle between atoms
    'atom'-'bond'-'angle' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::angleDelta(const AtomID &atom, const AtomID &bond,
                                const AtomID &angle) const
{
    return zmat.angleDelta(atom, bond, angle);
}

/** Return the maximum amount that the dihedral between atoms
    'atom'-'bond'-'angle'-'dihedral' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::dihedralDelta(const AtomID &atom, const AtomID &bond,
                                   const AtomID &angle, const AtomID &dihedral) const
{
    return zmat.dihedralDelta(atom, bond, angle, dihedral);
}

/** Return the maximum amount that the bond 'bond' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Length ZMatrixCoords::delta(const BondID &bond) const
{
    return this->bondDelta( bond.atom0(), bond.atom1() );
}

/** Return the maximum amount that the angle 'angle' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::delta(const AngleID &angle) const
{
    return this->angleDelta( angle.atom0(), angle.atom1(), angle.atom2() );
}

/** Return the maximum amount that the dihedral 'dihedral' should be changed

    \throw SireMol::missing_atom
    \throw SireMol::duplicate_atom
    \throw SireError::invalid_index
    \throw SireMove::zmatrix_error
*/
Angle ZMatrixCoords::delta(const DihedralID &dihedral) const
{
    return this->dihedralDelta( dihedral.atom0(), dihedral.atom1(),
                                dihedral.atom2(), dihedral.atom3() );
}

/** Return a z-matrix that only contains lines that involve the atoms
    that are in 'selection' */
ZMatrixCoords ZMatrixCoords::matchToSelection(const AtomSelection &selection) const
{
    this->rebuildCartesian();

    ZMatrixCoords new_zmat(*this);
    new_zmat.zmat = zmat.matchToSelection(selection);
    new_zmat.rebuildInternals();

    return new_zmat;
}

const char* ZMatrixCoords::typeName()
{
    return QMetaType::typeName( qMetaTypeId<ZMatrixCoords>() );
}
