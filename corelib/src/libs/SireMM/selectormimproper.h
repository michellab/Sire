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

#ifndef SIREMM_SELECTORMIMPROPER_H
#define SIREMM_SELECTORMIMPROPER_H

#include "SireMol/selectormol.h"
#include "SireMol/selectorm.hpp"

#include "selectorimproper.h"

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
class SelectorMImproper;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorMImproper&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorMImproper&);

namespace SireMM
{

/** Multi-molecule selector for impropers */
class SIREMM_EXPORT SelectorMImproper
    : public SireBase::ConcreteProperty<SelectorMImproper, SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorMImproper&);
friend QDataStream& ::operator>>(QDataStream&, SelectorMImproper&);

public:
    typedef QList<SelectorImproper>::const_iterator const_iterator;
    typedef QList<SelectorImproper>::const_iterator iterator;

    SelectorMImproper();
    SelectorMImproper(const Improper &view);

    SelectorMImproper(const SireMol::Molecules &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMImproper(const SireMol::MoleculeGroup &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());
    SelectorMImproper(const SireMol::MolGroupsBase &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMImproper(const SelectorImproper &Impropers);

    SelectorMImproper(const SireMol::SelectorMol &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    template<class U>
    SelectorMImproper(const SireMol::SelectorM<U> &other,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMImproper(const SelectorMImproper &impropers, const SireBase::Slice &slice);
    SelectorMImproper(const SelectorMImproper &impropers, const QList<qint64> &idxs);

    SelectorMImproper(const SireMol::SelectResult &mols,
                      const SireMol::ImproperID &Improper,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMImproper(const SireMol::SelectResult &mols,
                      const SireBase::PropertyMap &map=SireBase::PropertyMap());

    SelectorMImproper(const SireMol::SelectorM<SireMol::Atom> &atoms,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMImproper(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMImproper(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireMol::SelectorM<SireMol::Atom> &atoms2,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMImproper(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                      const SireMol::SelectorM<SireMol::Atom> &atoms1,
                      const SireMol::SelectorM<SireMol::Atom> &atoms2,
                      const SireMol::SelectorM<SireMol::Atom> &atoms3,
                      const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMImproper(const SelectorMImproper &other);

    virtual ~SelectorMImproper();

    static const char* typeName();

    virtual SelectorMImproper* clone() const
    {
        return new SelectorMImproper(*this);
    }

    SelectorMImproper& operator=(const SelectorMImproper &other);

    bool operator==(const SelectorMImproper &other) const;
    bool operator!=(const SelectorMImproper &other) const;

    Improper operator[](int i) const;
    SelectorMImproper operator[](const SireBase::Slice &slice) const;
    SelectorMImproper operator[](const QList<qint64> &idxs) const;
    SelectorMImproper operator[](const SireMol::ImproperID &id) const;

    Improper operator()(int i) const;
    SelectorMImproper operator()(const SireBase::Slice &slice) const;
    SelectorMImproper operator()(const QList<qint64> &idxs) const;
    SelectorMImproper operator()(const SireMol::ImproperID &id) const;

    bool isSelector() const;

    QList<SireMol::MolViewPtr> toList() const;
    SireMol::Molecules toMolecules() const;

    int count() const;
    int size() const;

    void update(const SireMol::Molecules &molecules);

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

    SelectorMImproper add(const SelectorMImproper &other) const;

    SelectorMImproper intersection(const SelectorMImproper &other) const;

    SelectorMImproper invert(const SireBase::PropertyMap &map) const;
    SelectorMImproper invert() const;

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

    QList<SireUnits::Dimension::GeneralUnit> energies() const;
    QList<SireUnits::Dimension::GeneralUnit> energies(
                            const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::GeneralUnit energy() const;
    SireUnits::Dimension::GeneralUnit energy(
                            const SireBase::PropertyMap &map) const;

    SireMol::SelectResult search(const QString &search_string) const;

    QList<SireMol::ImproperID> IDs() const;

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
    void _append(const Improper &improper);
    void _append(const SelectorImproper &impropers);

    /** The actual impropers */
    QList< SelectorImproper > imps;
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
SelectorMImproper::SelectorMImproper(const SireMol::SelectorM<T> &other,
                                     const SireBase::PropertyMap &map)
                  : SireBase::ConcreteProperty<SelectorMImproper,SireBase::Property>()
{
    if (not other.isEmpty())
    {
        this->imps.reserve(other.nMolecules());

        for (int i=0; i<other.nMolecules(); ++i)
        {
            this->imps.append(SelectorImproper(other.molecule(i), map));
        }
    }
}

} // end of namespace SireMol

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMM::SelectorMImproper)

SIRE_EXPOSE_CLASS(SireMM::SelectorMImproper)

SIRE_END_HEADER

#endif
