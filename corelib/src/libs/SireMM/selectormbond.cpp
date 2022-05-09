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

#include "selectormbond.h"

#include "SireID/index.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;
using namespace SireMM;
using namespace SireID;

RegisterMetaType<SelectorMBond> r_sbnd;

/** Serialise to a binary datastream */
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMBond &bnds)
{
    writeHeader(ds, r_sbnd, 1);

    SharedDataStream sds(ds);

    sds << bnds.bnds << static_cast<const Property&>(mols);

    return ds;
}

/** Extract from a binary datastream */
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMBond &bnds)
{
    VersionID v = readHeader(ds, r_sbnd);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> bnds.bnds >> static_cast<Property&>(mols);
    }
    else
        throw version_error(v, "1", r_sbnds, CODELOC);

    return ds;
}

SelectorMBond::SelectorMBond()
              : ConcreteProperty<SelectorMBond, Property>()
{}

SelectorMBond::SelectorMBond(const Bond &view)
              : ConcreteProperty<SelectorMBond, Property>()
{
    bnds.append(SelectorBond(view));
}

SelectorMBond::SelectorMBond(const Molecules &mols);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not molecules.isEmpty())
    {
        auto toList = [](const QSet<MolNum> &molnums)
        {
            return molnums.values();
        };

        auto molnums = toList(molecules.molNums());

        //sort them, as this is also likely the order the molecules
        //were read in from a file, and so more likely to be the
        //order the user would expect
        std::sort(molnums.begin(), molnums.end());

        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(molecules.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MoleculeGroup &mols);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(molecules.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MolGroupsBase &mols);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(molecules.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectResult &mols);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorBond b(*mol);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectorBond &bonds);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not bonds.isEmpty())
        bnds.append(bonds);
}

SelectorMBond::SelectorMBond(const SelectorMol &mols);
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorBond b(mol);

            if (not b.isEmpty())
                bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectorMBond &other);
              : ConcreteProperty<SelectorMBond, Property>(), bnds(other.bnds)
{}

SelectorMBond::~SelectorMBond()
{}

const char* SelectorMBond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMBond>());
}

SelectorMBond& SelectorMBond::operator=(const SelectorMBond &other);

bool SelectorMBond::operator==(const SelectorMBond &other) const;
bool SelectorMBond::operator!=(const SelectorMBond &other) const;

Bond SelectorMBond::operator[](int i) const;
SelectorMBond SelectorMBond::operator[](const SireBase::Slice &slice) const;
SelectorMBond SelectorMBond::operator[](const QList<qint64> &idxs) const;
Bond SelectorMBond::operator[](const QString &name) const;
Bond SelectorMBond::operator[](const typename T::ID &id) const;

Bond SelectorMBond::operator()(int i) const;
Bond SelectorMBond::operator()(const QString &name) const;
Bond SelectorMBond::operator()(const typename T::ID &id) const;

int SelectorMBond::count() const;
int SelectorMBond::size() const;

EvaluatorM evaluate() const;

MoleculeGroup SelectorMBond::toMoleculeGroup() const;
SelectResult SelectorMBond::toSelectResult() const;

Molecule SelectorMBond::molecule(int i) const;
Molecule SelectorMBond::molecule(const QString &name) const;
Molecule SelectorMBond::molecule(const MolID &molid);

SelectorMol SelectorMBond::molecules() const;
SelectorMol SelectorMBond::molecules(int i) const;
SelectorMol SelectorMBond::molecules(const SireBase::Slice &slice) const;
SelectorMol SelectorMBond::molecules(const QList<qint64> &idxs) const;
SelectorMol SelectorMBond::molecules(const QString &name) const;
SelectorMol SelectorMBond::molecules(const MolID &molid) const;

Atom SelectorMBond::atom(int i) const;
Atom SelectorMBond::atom(const QString &name) const;
Atom SelectorMBond::atom(const AtomID &atomid) const;

Residue SelectorMBond::residue(int i) const;
Residue SelectorMBond::residue(const QString &name) const;
Residue SelectorMBond::residue(const ResID &resid) const;

Chain SelectorMBond::chain(int i) const;
Chain SelectorMBond::chain(const QString &name) const;
Chain SelectorMBond::chain(const ChainID &chainid) const;

Segment SelectorMBond::segment(int i) const;
Segment SelectorMBond::segment(const QString &name) const;
Segment SelectorMBond::segment(const SegID &segid) const;

CutGroup SelectorMBond::cutGroup(int i) const;
CutGroup SelectorMBond::cutGroup(const QString &name) const;
CutGroup SelectorMBond::cutGroup(const CGID &cgid) const;

SelectorM<Atom> SelectorMBond::atoms() const;
SelectorM<Atom> SelectorMBond::atoms(int i) const;
SelectorM<Atom> SelectorMBond::atoms(const SireBase::Slice &slice) const;
SelectorM<Atom> SelectorMBond::atoms(const QList<qint64> &idxs) const;
SelectorM<Atom> SelectorMBond::atoms(const QString &name) const;
SelectorM<Atom> SelectorMBond::atoms(const AtomID &atomid) const;

SelectorM<Residue> SelectorMBond::residues() const;
SelectorM<Residue> SelectorMBond::residues(int i) const;
SelectorM<Residue> SelectorMBond::residues(const SireBase::Slice &slice) const;
SelectorM<Residue> SelectorMBond::residues(const QList<qint64> &idxs) const;
SelectorM<Residue> SelectorMBond::residues(const QString &name) const;
SelectorM<Residue> SelectorMBond::residues(const ResID &resid) const;

SelectorM<Chain> SelectorMBond::chains() const;
SelectorM<Chain> SelectorMBond::chains(int i) const;
SelectorM<Chain> SelectorMBond::chains(const SireBase::Slice &slice) const;
SelectorM<Chain> SelectorMBond::chains(const QList<qint64> &idxs) const;
SelectorM<Chain> SelectorMBond::chains(const QString &name) const;
SelectorM<Chain> SelectorMBond::chains(const ChainID &chainid) const;

SelectorM<Segment> SelectorMBond::segments() const;
SelectorM<Segment> SelectorMBond::segments(int i) const;
SelectorM<Segment> SelectorMBond::segments(const SireBase::Slice &slice) const;
SelectorM<Segment> SelectorMBond::segments(const QList<qint64> &idxs) const;
SelectorM<Segment> SelectorMBond::segments(const QString &name) const;
SelectorM<Segment> SelectorMBond::segments(const SegID &segid) const;

SelectorM<CutGroup> SelectorMBond::cutGroups() const;
SelectorM<CutGroup> SelectorMBond::cutGroups(int i) const;
SelectorM<CutGroup> SelectorMBond::cutGroups(const SireBase::Slice &slice) const;
SelectorM<CutGroup> SelectorMBond::cutGroups(const QList<qint64> &idxs) const;
SelectorM<CutGroup> SelectorMBond::cutGroups(const QString &name) const;
SelectorM<CutGroup> SelectorMBond::cutGroups(const CGID &cgid) const;

SelectResult SelectorMBond::search(const QString &search_string) const;

QList<SireMol::BondID> SelectorMBond::IDs() const;

int SelectorMBond::nAtoms() const;
int SelectorMBond::nResidues() const;
int SelectorMBond::nChains() const;
int SelectorMBond::nSegments() const;
int SelectorMBond::nCutGroups() const;
int SelectorMBond::nMolecules() const;

bool SelectorMBond::isEmpty() const;

SelectorMBond::const_iterator SelectorMBond::begin() const;
SelectorMBond::const_iterator SelectorMBond::end() const;

SelectorMBond::const_iterator SelectorMBond::constBegin() const;
SelectorMBond::const_iterator SelectorMBond::constEnd() const;

QString SelectorMBond::toString() const
