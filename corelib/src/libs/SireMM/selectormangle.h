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

#ifndef SIREMM_SELECTORMANGLE_H
#define SIREMM_SELECTORMANGLE_H

#include "SireMol/selectormol.h"
#include "SireMol/selectorm.hpp"

#include "selectorangle.h"

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
class SelectorMAngle;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::SelectorMAngle&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::SelectorMAngle&);

namespace SireMM
{

/** Multi-molecule selector for angles */
class SIREMM_EXPORT SelectorMAngle
    : public SireBase::ConcreteProperty<SelectorMAngle, SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const SelectorMAngle&);
friend QDataStream& ::operator>>(QDataStream&, SelectorMAngle&);

public:
    typedef QList<SelectorAngle>::const_iterator const_iterator;
    typedef QList<SelectorAngle>::const_iterator iterator;

    SelectorMAngle();
    SelectorMAngle(const Angle &view);

    SelectorMAngle(const SireMol::Molecules &mols,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorMAngle(const SireMol::MoleculeGroup &mols,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorMAngle(const SireMol::MolGroupsBase &mols,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorMAngle(const SireMol::SelectResult &mols,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());
    SelectorMAngle(const SelectorAngle &angles);

    SelectorMAngle(const SireMol::SelectorMol &mols,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    template<class U>
    SelectorMAngle(const SireMol::SelectorM<U> &other,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMAngle(const SelectorMAngle &angles, const SireBase::Slice &slice);
    SelectorMAngle(const SelectorMAngle &angles, const QList<qint64> &idxs);

    SelectorMAngle(const SireMol::SelectResult &mols,
                   const SireMol::AngleID &angle,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMAngle(const SireMol::SelectorM<SireMol::Atom> &atoms,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMAngle(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                   const SireMol::SelectorM<SireMol::Atom> &atoms1,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMAngle(const SireMol::SelectorM<SireMol::Atom> &atoms0,
                   const SireMol::SelectorM<SireMol::Atom> &atoms1,
                   const SireMol::SelectorM<SireMol::Atom> &atoms2,
                   const SireBase::PropertyMap &map = SireBase::PropertyMap());

    SelectorMAngle(const SelectorMAngle &other);

    virtual ~SelectorMAngle();

    static const char* typeName();

    virtual SelectorMAngle* clone() const
    {
        return new SelectorMAngle(*this);
    }

    SelectorMAngle& operator=(const SelectorMAngle &other);

    bool operator==(const SelectorMAngle &other) const;
    bool operator!=(const SelectorMAngle &other) const;

    Angle operator[](int i) const;
    SelectorMAngle operator[](const SireBase::Slice &slice) const;
    SelectorMAngle operator[](const QList<qint64> &idxs) const;
    SelectorMAngle operator[](const SireMol::AngleID &id) const;

    Angle operator()(int i) const;
    SelectorMAngle operator()(const SireBase::Slice &slice) const;
    SelectorMAngle operator()(const QList<qint64> &idxs) const;
    SelectorMAngle operator()(const SireMol::AngleID &id) const;

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

    SelectorMAngle add(const SelectorMAngle &other) const;

    SelectorMAngle intersection(const SelectorMAngle &other) const;

    SelectorMAngle invert(const SireBase::PropertyMap &map) const;
    SelectorMAngle invert() const;

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

    QList<SireMol::AngleID> IDs() const;

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
    void _append(const Angle &angle);
    void _append(const SelectorAngle &angles);

    /** The actual angles */
    QList< SelectorAngle > angs;
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
SelectorMAngle::SelectorMAngle(const SireMol::SelectorM<T> &other,
                               const SireBase::PropertyMap &map)
               : SireBase::ConcreteProperty<SelectorMAngle,SireBase::Property>()
{
    if (not other.isEmpty())
    {
        this->angs.reserve(other.nMolecules());

        for (int i=0; i<other.nMolecules(); ++i)
        {
            this->angs.append(SelectorAngle(other.molecule(i), map));
        }
    }
}

} // end of namespace SireMol

#endif // SIRE_SKIP_INLINE_FUNCTIONS

Q_DECLARE_METATYPE(SireMM::SelectorMAngle)

SIRE_EXPOSE_CLASS(SireMM::SelectorMAngle)

SIRE_END_HEADER

#endif
