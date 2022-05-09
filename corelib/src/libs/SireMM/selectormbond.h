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

#ifndef SIREMM_SELECTORMBOND_H
#define SIREMM_SELECTORMBOND_H

#include "selectorbond.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class SelectorMBond;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorMBond&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorMBond&);

namespace SireMM
{

/** Multi-molecule selector for bonds */
class SIREMM_EXPORT SelectorMBond
    : public SireBase::ConcreteProperty<SelectorMBond, SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorMBond&);
friend QDataStream& ::operator>>(QDataStream&, SelectorMBond&);

public:
    typedef QList<SelectorBond>::const_iterator const_iterator;
    typedef QList<SelectorBond>::const_iterator iterator;

    SelectorMBond();
    SelectorMBond(const Bond &view);
    SelectorMBond(const Molecules &mols);
    SelectorMBond(const MoleculeGroup &mols);
    SelectorMBond(const MolGroupsBase &mols);
    SelectorMBond(const SelectResult &mols);
    SelectorMBond(const SelectorBond &bonds);

    SelectorMBond(const SelectorMol &mols);

    template<class U>
    SelectorMBond(const SelectorM<U> &other);

    SelectorMBond(const SelectorMBond &other);

    virtual ~SelectorMBond();

    static const char* typeName();

    virtual SelectorMBond* clone() const
    {
        return new SelectorMBond(*this);
    }

    SelectorMBond& operator=(const SelectorMBond &other);

    bool operator==(const SelectorMBond &other) const;
    bool operator!=(const SelectorMBond &other) const;

    Bond operator[](int i) const;
    SelectorMBond operator[](const SireBase::Slice &slice) const;
    SelectorMBond operator[](const QList<qint64> &idxs) const;
    Bond operator[](const QString &name) const;
    Bond operator[](const typename T::ID &id) const;

    Bond operator()(int i) const;
    Bond operator()(const QString &name) const;
    Bond operator()(const typename T::ID &id) const;

    int count() const;
    int size() const;

    EvaluatorM evaluate() const;

    MoleculeGroup toMoleculeGroup() const;
    SelectResult toSelectResult() const;

    Molecule molecule(int i) const;
    Molecule molecule(const QString &name) const;
    Molecule molecule(const MolID &molid);

    SelectorMol molecules() const;
    SelectorMol molecules(int i) const;
    SelectorMol molecules(const SireBase::Slice &slice) const;
    SelectorMol molecules(const QList<qint64> &idxs) const;
    SelectorMol molecules(const QString &name) const;
    SelectorMol molecules(const MolID &molid) const;

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

    SelectResult search(const QString &search_string) const;

    QList<SireMol::BondID> IDs() const;

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
    void _append(const Bond &bond);

    /** The actual bonds */
    QList< SelectorBond > bnds;
};

} // end of namespace SireMM

#include "SireMol/selectorm.hpp"

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

namespace SireMM
{

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMBond::SelectorMBond(const SelectorM<T> &other)
              : SireBase::ConcreteProperty<SelectorMBond,SireBase::Property>()
{
    if (not other.isEmpty())
    {
        this->bnds.reserve(other.nMolecules());

        for (int i=0; i<other.nMolecules(); ++i)
        {
            this->bnds.append(SelectorBond(other.molecule(i)));
        }
    }
}

} // end of namespace SireMol

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMM::SelectorMBond)

SIRE_EXPOSE_CLASS(SireMM::SelectorMBond)

SIRE_END_HEADER

#endif
