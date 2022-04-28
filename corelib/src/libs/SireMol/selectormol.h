/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMOL_SELECTORMOL_H
#define SIREMOL_SELECTORMOL_H

#include "molecule.h"
#include "atom.h"
#include "residue.h"
#include "chain.h"
#include "segment.h"
#include "cutgroup.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class SelectorMol;

    template<class T>
    class SelectorM;

    class Molecules;
    class MoleculeGroup;
    class MolGroupsBase;

    class AtomID;
    class ResID;
    class ChainID;
    class SegID;
    class CGID;

    class MolID;
    class MolName;
    class MolNum;
    class MolIdx;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::SelectorMol&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::SelectorMol&);

namespace SireMol
{

/** This class provides a Selector<T>-type interface to a
    collection of molecules.
*/
class SIREMOL_EXPORT SelectorMol
    : public SireBase::ConcreteProperty<SelectorMol,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const SelectorMol&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, SelectorMol&);

public:
    typedef QList<Molecule>::const_iterator const_iterator;
    typedef QList<Molecule>::const_iterator iterator;

    SelectorMol();
    SelectorMol(const MoleculeView &mol);
    SelectorMol(const Molecules &mols);
    SelectorMol(const MoleculeGroup &mols);
    SelectorMol(const MolGroupsBase &mols);
    SelectorMol(const SelectorMol &other);

    virtual ~SelectorMol();

    static const char* typeName();

    virtual const char* what() const
    {
        return SelectorMol::typeName();
    }

    virtual SelectorMol* clone() const
    {
        return new SelectorMol(*this);
    }

    SelectorMol& operator=(const SelectorMol &other);

    bool operator==(const SelectorMol &other) const;
    bool operator!=(const SelectorMol &other) const;

    Molecule operator[](int i) const;
    Molecule operator[](const QString &name) const;
    Molecule operator[](const MolIdx &molidx) const;
    Molecule operator[](const MolName &molname) const;
    Molecule operator[](const MolNum &molnum) const;
    Molecule operator[](const MolID &molid) const;

    int count() const;
    int size() const;

    Molecule molecule(int i) const;
    Molecule molecule(const QString &name) const;
    Molecule molecule(const MolIdx &molidx) const;
    Molecule molecule(const MolName &molname) const;
    Molecule molecule(const MolNum &molnum) const;
    Molecule molecule(const MolID &molid);

    SelectorMol molecules() const;
    SelectorMol molecules(int i) const;
    SelectorMol molecules(const SireBase::Slice &slice) const;
    SelectorMol molecules(const QList<qint64> &idxs) const;
    SelectorMol molecules(const QString &name) const;
    SelectorMol molecules(const MolIdx &molidx) const;
    SelectorMol molecules(const MolName &molname) const;
    SelectorMol molecules(const MolNum &molnum) const;
    SelectorMol molecules(const MolID &molid) const;

    bool contains(const MolNum &molnum) const;

    Atom atom(int i) const;
    Atom atom(const QString &name) const;
    Atom atom(const AtomID &atomid) const;

    Residue residue(int i) const;
    Residue residue(const QString &name) const;
    Residue residue(const ResID &resid) const;

    Chain chain(int i) const;
    Chain chain(const QString &name) const;
    Chain chain(const ChainID &chainid) const;

    Segment segment(int i) const;
    Segment segment(const QString &name) const;
    Segment segment(const SegID &segid) const;

    CutGroup cutGroup(int i) const;
    CutGroup cutGroup(const QString &name) const;
    CutGroup cutGroup(const CGID &cgid) const;

    SelectorM<Atom> atoms() const;
    SelectorM<Atom> atoms(int i) const;
    SelectorM<Atom> atoms(const SireBase::Slice &slice) const;
    SelectorM<Atom> atoms(const QList<qint64> &idxs) const;
    SelectorM<Atom> atoms(const QString &name) const;
    SelectorM<Atom> atoms(const AtomID &atomid) const;

    SelectorM<Residue> residues() const;
    SelectorM<Residue> residues(int i) const;
    SelectorM<Residue> residues(const SireBase::Slice &slice) const;
    SelectorM<Residue> residues(const QList<qint64> &idxs) const;
    SelectorM<Residue> residues(const QString &name) const;
    SelectorM<Residue> residues(const ResID &resid) const;

    SelectorM<Chain> chains() const;
    SelectorM<Chain> chains(int i) const;
    SelectorM<Chain> chains(const SireBase::Slice &slice) const;
    SelectorM<Chain> chains(const QList<qint64> &idxs) const;
    SelectorM<Chain> chains(const QString &name) const;
    SelectorM<Chain> chains(const ChainID &chainid) const;

    SelectorM<Segment> segments() const;
    SelectorM<Segment> segments(int i) const;
    SelectorM<Segment> segments(const SireBase::Slice &slice) const;
    SelectorM<Segment> segments(const QList<qint64> &idxs) const;
    SelectorM<Segment> segments(const QString &name) const;
    SelectorM<Segment> segments(const SegID &segid) const;

    SelectorM<CutGroup> cutGroups() const;
    SelectorM<CutGroup> cutGroups(int i) const;
    SelectorM<CutGroup> cutGroups(const SireBase::Slice &slice) const;
    SelectorM<CutGroup> cutGroups(const QList<qint64> &idxs) const;
    SelectorM<CutGroup> cutGroups(const QString &name) const;
    SelectorM<CutGroup> cutGroups(const CGID &cgid) const;

    QVector<MolNum> molNums() const;

    QList<MolNum> IDs() const;
    QList<MolIdx> indexes() const;
    QList<MolNum> numbers() const;
    QList<MolName> names() const;

    int nAtoms() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;
    int nCutGroups() const;
    int nMolecules() const;

    bool isEmpty() const;

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;

    virtual QString toString() const;

protected:
    /** The actual list of molecules */
    QList<Molecule> mols;
};

}

Q_DECLARE_METATYPE(SireMol::SelectorMol);
Q_DECLARE_METATYPE(SireMol::SelectorM<SireMol::Atom>);
Q_DECLARE_METATYPE(SireMol::SelectorM<SireMol::Residue>);
Q_DECLARE_METATYPE(SireMol::SelectorM<SireMol::Chain>);
Q_DECLARE_METATYPE(SireMol::SelectorM<SireMol::Segment>);
Q_DECLARE_METATYPE(SireMol::SelectorM<SireMol::CutGroup>);

SIRE_EXPOSE_CLASS(SireMol::SelectorMol)

SIRE_EXPOSE_ALIAS( SireMol::SelectorM<SireMol::Atom>, SireMol::SelectorM_Atom_ )
SIRE_EXPOSE_ALIAS( SireMol::SelectorM<SireMol::Residue>, SireMol::SelectorM_Residue_ )
SIRE_EXPOSE_ALIAS( SireMol::SelectorM<SireMol::Chain>, SireMol::SelectorM_Chain_ )
SIRE_EXPOSE_ALIAS( SireMol::SelectorM<SireMol::Segment>, SireMol::SelectorM_Segment_ )
SIRE_EXPOSE_ALIAS( SireMol::SelectorM<SireMol::CutGroup>, SireMol::SelectorM_CutGroup_ )

// include this here so that SelectorM is always instantiatable if
// SelectorMol has been included
#include "selectorm.hpp"

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "selectorm.hpp"

template class SireMol::SelectorM<SireMol::Atom>;
template class SireMol::SelectorM<SireMol::Residue>;
template class SireMol::SelectorM<SireMol::Chain>;
template class SireMol::SelectorM<SireMol::Segment>;
template class SireMol::SelectorM<SireMol::CutGroup>;

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
