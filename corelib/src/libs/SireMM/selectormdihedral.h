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

#ifndef SIREMM_SELECTORMDIHEDRAL_H
#define SIREMM_SELECTORMDIHEDRAL_H

#include "SireMol/selectormol.h"
#include "SireMol/selectorm.hpp"

#include "selectordihedral.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Atom;
class Residue;
class Chain;
class Segment;
class CutGroup;
class Molecule;
class MoleculeGroup;
class Molecules;
class SelectResult;

class SelectorMol;
class EvaluatorM;

}

namespace SireMM
{
class SelectorMDihedral;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorMDihedral&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorMDihedral&);

namespace SireMM
{

/** Multi-molecule selector for dihedrals */
class SIREMM_EXPORT SelectorMDihedral
    : public SireBase::ConcreteProperty<SelectorMDihedral, SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorMDihedral&);
friend QDataStream& ::operator>>(QDataStream&, SelectorMDihedral&);

public:
    typedef QList<SelectorDihedral>::const_iterator const_iterator;
    typedef QList<SelectorDihedral>::const_iterator iterator;

    SelectorMDihedral();
    SelectorMDihedral(const Dihedral &view);

    SelectorMDihedral(const SireMol::Molecules &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMDihedral(const SireMol::MoleculeGroup &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMDihedral(const SireMol::MolGroupsBase &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMDihedral(const SireMol::SelectResult &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMDihedral(const SelectorDihedral &Dihedrals);

    SelectorMDihedral(const SireMol::SelectorMol &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    template<class U>
    SelectorMDihedral(const SireMol::SelectorM<U> &other,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMDihedral(const SelectorMDihedral &dihedrals, const SireBase::Slice &slice);
    SelectorMDihedral(const SelectorMDihedral &dihedrals, const QList<qint64> &idxs);

    SelectorMDihedral(const SireMol::SelectResult &mols,
                      const SireMol::DihedralID &Dihedral,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMDihedral(const SireMol::SelectorM<SireMol::Atom> &atoms,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMDihedral(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMDihedral(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireMol::SelectorM<SireMol::Atom> &atoms2,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMDihedral(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireMol::SelectorM<SireMol::Atom> &atoms2,
                      const SireMol::SelectorM<SireMol::Atom> &atoms3,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMDihedral(const SelectorMDihedral &other);

    virtual ~SelectorMDihedral();

    static const char* typeName();

    virtual SelectorMDihedral* clone() const
    {
        return new SelectorMDihedral(*this);
    }

    SelectorMDihedral& operator=(const SelectorMDihedral &other);

    bool operator==(const SelectorMDihedral &other) const;
    bool operator!=(const SelectorMDihedral &other) const;

    Dihedral operator[](int i) const;
    SelectorMDihedral operator[](const SireBase::Slice &slice) const;
    SelectorMDihedral operator[](const QList<qint64> &idxs) const;
    SelectorMDihedral operator[](const SireMol::DihedralID &id) const;

    Dihedral operator()(int i) const;
    SelectorMDihedral operator()(const SireBase::Slice &slice) const;
    SelectorMDihedral operator()(const QList<qint64> &idxs) const;
    SelectorMDihedral operator()(const SireMol::DihedralID &id) const;

    QList<SireMol::MolViewPtr> toList() const;

    int count() const;
    int size() const;

    SireMol::EvaluatorM evaluate() const;

    SireMol::MoleculeGroup toMoleculeGroup() const;
    SireMol::SelectResult toSelectResult() const;

    SireMol::Molecule molecule(int i) const;
    SireMol::Molecule molecule(const QString &name) const;
    SireMol::Molecule molecule(const SireMol::MolID &molid);

    SireMol::SelectorMol molecules() const;
    SireMol::SelectorMol molecules(int i) const;
    SireMol::SelectorMol molecules(const SireBase::Slice &slice) const;
    SireMol::SelectorMol molecules(const QList<qint64> &idxs) const;
    SireMol::SelectorMol molecules(const QString &name) const;
    SireMol::SelectorMol molecules(const SireMol::MolID &molid) const;

    SireMol::Atom atom(int i) const;
    SireMol::Atom atom(const QString &name) const;
    SireMol::Atom atom(const SireMol::AtomID &atomid) const;

    SireMol::Residue residue(int i) const;
    SireMol::Residue residue(const QString &name) const;
    SireMol::Residue residue(const SireMol::ResID &resid) const;

    SireMol::Chain chain(int i) const;
    SireMol::Chain chain(const QString &name) const;
    SireMol::Chain chain(const SireMol::ChainID &chainid) const;

    SireMol::Segment segment(int i) const;
    SireMol::Segment segment(const QString &name) const;
    SireMol::Segment segment(const SireMol::SegID &segid) const;

    SireMol::CutGroup cutGroup(int i) const;
    SireMol::CutGroup cutGroup(const QString &name) const;
    SireMol::CutGroup cutGroup(const SireMol::CGID &cgid) const;

    SireMol::SelectorM<SireMol::Atom> atoms() const;
    SireMol::SelectorM<SireMol::Atom> atoms(int i) const;
    SireMol::SelectorM<SireMol::Atom> atoms(const SireBase::Slice &slice) const;
    SireMol::SelectorM<SireMol::Atom> atoms(const QList<qint64> &idxs) const;
    SireMol::SelectorM<SireMol::Atom> atoms(const QString &name) const;
    SireMol::SelectorM<SireMol::Atom> atoms(const SireMol::AtomID &atomid) const;

    SireMol::SelectorM<SireMol::Residue> residues() const;
    SireMol::SelectorM<SireMol::Residue> residues(int i) const;
    SireMol::SelectorM<SireMol::Residue> residues(const SireBase::Slice &slice) const;
    SireMol::SelectorM<SireMol::Residue> residues(const QList<qint64> &idxs) const;
    SireMol::SelectorM<SireMol::Residue> residues(const QString &name) const;
    SireMol::SelectorM<SireMol::Residue> residues(const SireMol::ResID &resid) const;

    SireMol::SelectorM<SireMol::Chain> chains() const;
    SireMol::SelectorM<SireMol::Chain> chains(int i) const;
    SireMol::SelectorM<SireMol::Chain> chains(const SireBase::Slice &slice) const;
    SireMol::SelectorM<SireMol::Chain> chains(const QList<qint64> &idxs) const;
    SireMol::SelectorM<SireMol::Chain> chains(const QString &name) const;
    SireMol::SelectorM<SireMol::Chain> chains(const SireMol::ChainID &chainid) const;

    SireMol::SelectorM<SireMol::Segment> segments() const;
    SireMol::SelectorM<SireMol::Segment> segments(int i) const;
    SireMol::SelectorM<SireMol::Segment> segments(const SireBase::Slice &slice) const;
    SireMol::SelectorM<SireMol::Segment> segments(const QList<qint64> &idxs) const;
    SireMol::SelectorM<SireMol::Segment> segments(const QString &name) const;
    SireMol::SelectorM<SireMol::Segment> segments(const SireMol::SegID &segid) const;

    SireMol::SelectorM<SireMol::CutGroup> cutGroups() const;
    SireMol::SelectorM<SireMol::CutGroup> cutGroups(int i) const;
    SireMol::SelectorM<SireMol::CutGroup> cutGroups(const SireBase::Slice &slice) const;
    SireMol::SelectorM<SireMol::CutGroup> cutGroups(const QList<qint64> &idxs) const;
    SireMol::SelectorM<SireMol::CutGroup> cutGroups(const QString &name) const;
    SireMol::SelectorM<SireMol::CutGroup> cutGroups(const SireMol::CGID &cgid) const;

    SelectorMDihedral add(const SelectorMDihedral &other) const;

    SelectorMDihedral intersection(const SelectorMDihedral &other) const;

    SelectorMDihedral invert(const SireBase::PropertyMap &map) const;
    SelectorMDihedral invert() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    QList<SireBase::Properties> properties() const;

    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key) const;
    QList<SireBase::PropertyPtr> property(const SireBase::PropertyName &key,
                                          const Property &default_value) const;

    QList<SireUnits::Dimension::Angle> sizes() const;
    QList<SireUnits::Dimension::Angle> sizes(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::Angle> measures() const;
    QList<SireUnits::Dimension::Angle> measures(const SireBase::PropertyMap &map) const;

    QList<SireCAS::Expression> potentials() const;
    QList<SireCAS::Expression> potentials(const SireBase::PropertyMap &map) const;

    QList<SireUnits::Dimension::MolarEnergy> energies() const;
    QList<SireUnits::Dimension::MolarEnergy> energies(
                            const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

    SireMol::SelectResult search(const QString &search_string) const;

    QList<SireMol::DihedralID> IDs() const;

    int nAtoms() const;
    int nResidues() const;
    int nChains() const;
    int nSegments() const;
    int nCutGroups() const;
    int nMolecules() const;

    bool isEmpty() const;

    int nFrames() const;
    int nFrames(const SireBase::PropertyMap &map) const;

    void loadFrame(int frame);
    void saveFrame(int frame);
    void saveFrame();
    void deleteFrame(int frame);

    void loadFrame(int frame, const SireBase::PropertyMap &map);
    void saveFrame(int frame, const SireBase::PropertyMap &map);
    void saveFrame(const SireBase::PropertyMap &map);
    void deleteFrame(int frame, const SireBase::PropertyMap &map);

    const_iterator begin() const;
    const_iterator end() const;

    const_iterator constBegin() const;
    const_iterator constEnd() const;

    virtual QString toString() const;

protected:
    void _append(const Dihedral &dihedral);
    void _append(const SelectorDihedral &dihedrals);

    /** The actual dihedrals */
    QList< SelectorDihedral > dihs;
};

} // end of namespace SireMM

#include "SireMol/selectorm.hpp"
#include "SireMol/atom.h"
#include "SireMol/residue.h"
#include "SireMol/chain.h"
#include "SireMol/segment.h"
#include "SireMol/cutgroup.h"

#ifndef SIRE_SKIP_INLINE_FUNCTIONS

namespace SireMM
{

template<class T>
SIRE_OUTOFLINE_TEMPLATE
SelectorMDihedral::SelectorMDihedral(const SireMol::SelectorM<T> &other,
                                     const SireBase::PropertyMap &map)
                  : SireBase::ConcreteProperty<SelectorMDihedral,SireBase::Property>()
{
    if (not other.isEmpty())
    {
        this->dihs.reserve(other.nMolecules());

        for (int i=0; i<other.nMolecules(); ++i)
        {
            this->dihs.append(SelectorDihedral(other.molecule(i), map));
        }
    }
}

} // end of namespace SireMol

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMM::SelectorMDihedral)

SIRE_EXPOSE_CLASS(SireMM::SelectorMDihedral)

SIRE_END_HEADER

#endif
