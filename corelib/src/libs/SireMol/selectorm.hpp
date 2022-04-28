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

#ifndef SIREMOL_SELECTORM_HPP
#define SIREMOL_SELECTORM_HPP

#include "selectormol.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
template<class T>
class SelectorM;

class SelectorMol;
}

template<class T>
QDataStream& operator<<(QDataStream&, const SelectorM<T>&);

template<class T>
QDataStream& operator>>(QDataStream&, SelectorM<T>&);

namespace SireMol
{

/** This is an analogue of the Selector<T> class that is designed
    to hold views from multiple molecules
*/
template<class T>
class SIREMOL_EXPORT SelectorM
    : public SireBase::ConcreteProperty<SelectorM<T>,SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<<>(QDataStream&, const SelectorM<T>&);
friend SIREMOL_EXPORT QDataStream& ::operator>><>(QDataStream&, SelectorM<T>&);

public:
    SelectorM();
    SelectorM(const T &view);

    SelectorM(const SelectorMol &mols);
    SelectorM(const SelectorMol &mols, const SireBase::Slice &slice);
    SelectorM(const SelectorMol &mols, const QList<qint64> &idxs);
    SelectorM(const SelectorMol &mols, const QString &name);
    SelectorM(const SelectorMol &mols, const T::ID &id);

    template<class U>
    SelectorM(const SelectorM<U> &other);
    template<class U>
    SelectorM(const SelectorM<U> &other, const SireBase::Slice &slice);
    template<class U>
    SelectorM(const SelectorM<U> &other, const QList<qint64> &idxs);
    template<class U>
    SelectorM(const SelectorM<U> &other, const QString &name);
    template<class U>
    SelectorM(const SelectorM<U> &other, const T::ID &id);

    SelectorM(const SelectorM<T> &other);

    virtual ~SelectorM();

    static const char* typeName();

    virtual SelectorM<T>* clone() const
    {
        return new SelectorM<T>(*this);
    }

    SelectorM<T>& operator=(const SelectorM<T> &other);

    bool operator==(const SelectorM<T> &other) const;
    bool operator!=(const SelectorM<T> &other) const;

    T operator[](int i) const;
    T operator[](const QString &name) const;
    T operator[](const T::ID &id) const;

    int count() const;
    int size() const;

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

    QList<T::ID> IDs() const;
    QList<T::Index> indexes() const;
    QList<T::Number> numbers() const;
    QList<T::Name> names() const;

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
};

}

SIRE_END_HEADER

#endif
