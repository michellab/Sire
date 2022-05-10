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

    sds << bnds.bnds << static_cast<const Property&>(bnds);

    return ds;
}

/** Extract from a binary datastream */
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMBond &bnds)
{
    VersionID v = readHeader(ds, r_sbnd);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> bnds.bnds >> static_cast<Property&>(bnds);
    }
    else
        throw version_error(v, "1", r_sbnd, CODELOC);

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

SelectorMBond::SelectorMBond(const Molecules &mols)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        auto toList = [](const QSet<MolNum> &molnums)
        {
            return molnums.values();
        };

        auto molnums = toList(mols.molNums());

        //sort them, as this is also likely the order the molecules
        //were read in from a file, and so more likely to be the
        //order the user would expect
        std::sort(molnums.begin(), molnums.end());

        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(mols.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MoleculeGroup &mols)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(mols.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MolGroupsBase &mols)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(mols.at(molnum));

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectResult &mols)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorBond b;

            if (mol->isA<SelectorBond>())
                b = mol->asA<SelectorBond>();
            else
                b = SelectorBond(*mol);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectResult &mols, const BondID &bond)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            try
            {
                auto b = SelectorBond(*mol, bond);

                if (not b.isEmpty())
                    this->bnds.append(b);
            }
            catch(...)
            {}
        }
    }
}

SelectorMBond::SelectorMBond(const SelectorBond &bonds)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not bonds.isEmpty())
        bnds.append(bonds);
}

SelectorMBond::SelectorMBond(const SelectorMol &mols)
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

void SelectorMBond::_append(const Bond &bond)
{
    if (this->bnds.isEmpty())
    {
        this->bnds.append(SelectorBond(bond));
    }
    else if (this->bnds.last().data().number() != bond.data().number())
    {
        // new molecule
        this->bnds.append(SelectorBond(bond));
    }
    else
    {
        // a new view in the current molecule
        this->bnds.last() = this->bnds.last().add(bond);
    }
}

SelectorMBond::SelectorMBond(const SelectorMBond &bonds,
                             const SireBase::Slice &slice)
              : SireBase::ConcreteProperty<SelectorMBond,Property>()
{
    for (auto it = slice.begin(bonds.count());
         not it.atEnd(); it.next())
    {
        this->_append(bonds[it.value()]);
    }
}

SelectorMBond::SelectorMBond(const SelectorMBond &bonds,
                             const QList<qint64> &idxs)
              : SireBase::ConcreteProperty<SelectorMBond,Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(bonds[idx]);
    }
}

SelectorMBond::SelectorMBond(const SelectorMBond &other)
              : ConcreteProperty<SelectorMBond, Property>(), bnds(other.bnds)
{}

SelectorMBond::~SelectorMBond()
{}

const char* SelectorMBond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMBond>());
}

SelectorMBond& SelectorMBond::operator=(const SelectorMBond &other)
{
    if (this != &other)
    {
        bnds = other.bnds;
        Property::operator=(other);
    }

    return *this;
}

bool SelectorMBond::operator==(const SelectorMBond &other) const
{
    return bnds == other.bnds;
}

bool SelectorMBond::operator!=(const SelectorMBond &other) const
{
    return not operator==(other);
}

Bond SelectorMBond::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    for (const auto &b : bnds)
    {
        if (i < b.count())
        {
            return b(i);
        }
        else
        {
            i -= b.count();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Bond();
}

SelectorMBond SelectorMBond::operator[](const SireBase::Slice &slice) const
{
    return SelectorMBond(*this, slice);
}

SelectorMBond SelectorMBond::operator[](const QList<qint64> &idxs) const
{
    return SelectorMBond(*this, idxs);
}

SelectorMBond SelectorMBond::operator[](const BondID &id) const
{
    SelectorMBond ret;

    for (const auto &b : bnds)
    {
        try
        {
            auto r = b(id);

            if (not r.isEmpty())
            {
                ret.bnds.append(r);
            }
        }
        catch(...)
        {}
    }

    return ret;
}

Bond SelectorMBond::operator()(int i) const
{
    return this->operator[](i);
}

SelectorMBond SelectorMBond::operator()(const SireBase::Slice &slice) const
{
    return this->operator[](slice);
}

SelectorMBond SelectorMBond::operator()(const QList<qint64> &idxs) const
{
    return this->operator[](idxs);
}

SelectorMBond SelectorMBond::operator()(const BondID &id) const
{
    return this->operator[](id);
}

int SelectorMBond::count() const
{
    int n = 0;

    for (const auto &b : bnds)
    {
        n += b.count();
    }

    return n;
}

int SelectorMBond::size() const
{
    return this->count();
}

EvaluatorM SelectorMBond::evaluate() const
{
    return EvaluatorM(this->atoms());
}

MoleculeGroup SelectorMBond::toMoleculeGroup() const
{
    MoleculeGroup grp;

    for (const auto &b : this->bnds)
    {
        grp.add(b);
    }

    return grp;

}

SelectResult SelectorMBond::toSelectResult() const
{
    QList<MolViewPtr> r;

    for (const auto &b : bnds)
    {
        r.append(b);
    }

    return SelectResult(r);
}

Molecule SelectorMBond::molecule(int i) const
{
    return this->molecules().molecule(i);
}

Molecule SelectorMBond::molecule(const QString &name) const
{
    return this->molecules().molecule(name);
}

Molecule SelectorMBond::molecule(const MolID &molid)
{
    return this->molecules().molecule(molid);
}

SelectorMol SelectorMBond::molecules() const
{
    QList<Molecule> mols;

    for (const auto &b : this->bnds)
    {
        mols.append(b.molecule());
    }

    return SelectorMol(mols);
}

SelectorMol SelectorMBond::molecules(int i) const
{
    return this->molecules().molecules(i);
}

SelectorMol SelectorMBond::molecules(const SireBase::Slice &slice) const
{
    return this->molecules().molecules(slice);
}

SelectorMol SelectorMBond::molecules(const QList<qint64> &idxs) const
{
    return this->molecules().molecules(idxs);
}

SelectorMol SelectorMBond::molecules(const QString &name) const
{
    return this->molecules().molecules(name);
}

SelectorMol SelectorMBond::molecules(const MolID &molid) const
{
    return this->molecules().molecules(molid);
}

Atom SelectorMBond::atom(int i) const
{
    return this->atoms()(i);
}

Atom SelectorMBond::atom(const QString &name) const
{
    return this->atoms()(name);
}

Atom SelectorMBond::atom(const AtomID &atomid) const
{
    return this->atoms()(atomid);
}

Residue SelectorMBond::residue(int i) const
{
    return this->residues()(i);
}

Residue SelectorMBond::residue(const QString &name) const
{
    return this->residues()(name);
}

Residue SelectorMBond::residue(const ResID &resid) const
{
    return this->residues()(resid);
}

Chain SelectorMBond::chain(int i) const
{
    return this->chains()(i);
}

Chain SelectorMBond::chain(const QString &name) const
{
    return this->chains()(name);
}

Chain SelectorMBond::chain(const ChainID &chainid) const
{
    return this->chains()(chainid);
}

Segment SelectorMBond::segment(int i) const
{
    return this->segments()(i);
}

Segment SelectorMBond::segment(const QString &name) const
{
    return this->segments()(name);
}

Segment SelectorMBond::segment(const SegID &segid) const
{
    return this->segments()(segid);
}

CutGroup SelectorMBond::cutGroup(int i) const
{
    return this->cutGroups()(i);
}

CutGroup SelectorMBond::cutGroup(const QString &name) const
{
    return this->cutGroups()(name);
}

CutGroup SelectorMBond::cutGroup(const CGID &cgid) const
{
    return this->cutGroups()(cgid);
}

SelectorM<Atom> SelectorMBond::atoms() const
{
    QList< Selector<Atom> > ret;

    for (const auto &b : this->bnds)
    {
        ret.append(b.atoms());
    }

    return SelectorM<Atom>(ret);
}

SelectorM<Atom> SelectorMBond::atoms(int i) const
{
    return this->atoms().atoms(i);
}

SelectorM<Atom> SelectorMBond::atoms(const SireBase::Slice &slice) const
{
    return this->atoms().atoms(slice);
}

SelectorM<Atom> SelectorMBond::atoms(const QList<qint64> &idxs) const
{
    return this->atoms().atoms(idxs);
}

SelectorM<Atom> SelectorMBond::atoms(const QString &name) const
{
    return this->atoms().atoms(name);
}

SelectorM<Atom> SelectorMBond::atoms(const AtomID &atomid) const
{
    return this->atoms().atoms(atomid);
}

SelectorM<Residue> SelectorMBond::residues() const
{
    QList< Selector<Residue> > ret;

    for (const auto &b : this->bnds)
    {
        ret.append(b.residues());
    }

    return SelectorM<Residue>(ret);
}

SelectorM<Residue> SelectorMBond::residues(int i) const
{
    return this->residues().residues(i);
}

SelectorM<Residue> SelectorMBond::residues(const SireBase::Slice &slice) const
{
    return this->residues().residues(slice);
}

SelectorM<Residue> SelectorMBond::residues(const QList<qint64> &idxs) const
{
    return this->residues().residues(idxs);
}

SelectorM<Residue> SelectorMBond::residues(const QString &name) const
{
    return this->residues().residues(name);
}

SelectorM<Residue> SelectorMBond::residues(const ResID &resid) const
{
    return this->residues().residues(resid);
}

SelectorM<Chain> SelectorMBond::chains() const
{
    QList< Selector<Chain> > ret;

    for (const auto &b : this->bnds)
    {
        ret.append(b.chains());
    }

    return SelectorM<Chain>(ret);
}

SelectorM<Chain> SelectorMBond::chains(int i) const
{
    return this->chains().chains(i);
}

SelectorM<Chain> SelectorMBond::chains(const SireBase::Slice &slice) const
{
    return this->chains().chains(slice);
}

SelectorM<Chain> SelectorMBond::chains(const QList<qint64> &idxs) const
{
    return this->chains().chains(idxs);
}

SelectorM<Chain> SelectorMBond::chains(const QString &name) const
{
    return this->chains().chains(name);
}

SelectorM<Chain> SelectorMBond::chains(const ChainID &chainid) const
{
    return this->chains().chains(chainid);
}

SelectorM<Segment> SelectorMBond::segments() const
{
    QList< Selector<Segment> > ret;

    for (const auto &b : this->bnds)
    {
        ret.append(b.segments());
    }

    return SelectorM<Segment>(ret);
}

SelectorM<Segment> SelectorMBond::segments(int i) const
{
    return this->segments().segments(i);
}

SelectorM<Segment> SelectorMBond::segments(const SireBase::Slice &slice) const
{
    return this->segments().segments(slice);
}

SelectorM<Segment> SelectorMBond::segments(const QList<qint64> &idxs) const
{
    return this->segments().segments(idxs);
}

SelectorM<Segment> SelectorMBond::segments(const QString &name) const
{
    return this->segments().segments(name);
}

SelectorM<Segment> SelectorMBond::segments(const SegID &segid) const
{
    return this->segments().segments(segid);
}

SelectorM<CutGroup> SelectorMBond::cutGroups() const
{
    QList< Selector<CutGroup> > ret;

    for (const auto &b : this->bnds)
    {
        ret.append(b.cutGroups());
    }

    return SelectorM<CutGroup>(ret);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(int i) const
{
    return this->cutGroups().cutGroups(i);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const SireBase::Slice &slice) const
{
    return this->cutGroups().cutGroups(slice);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const QList<qint64> &idxs) const
{
    return this->cutGroups().cutGroups(idxs);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const QString &name) const
{
    return this->cutGroups().cutGroups(name);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const CGID &cgid) const
{
    return this->cutGroups().cutGroups(cgid);
}

SelectResult SelectorMBond::search(const QString &search_string) const
{
    return this->toSelectResult().search(search_string);
}

QList<BondID> SelectorMBond::IDs() const
{
    QList<BondID> ret;

    for (const auto &b : this->bnds)
    {
        ret += b.IDs();
    }

    return ret;
}

int SelectorMBond::nAtoms() const
{
    int n = 0;

    for (const auto &b : this->bnds)
    {
        n += b.nAtoms();
    }

    return n;
}

int SelectorMBond::nResidues() const
{
    int n = 0;

    for (const auto &b : this->bnds)
    {
        n += b.nResidues();
    }

    return n;
}

int SelectorMBond::nChains() const
{
    int n = 0;

    for (const auto &b : this->bnds)
    {
        n += b.nChains();
    }

    return n;
}

int SelectorMBond::nSegments() const
{
    int n = 0;

    for (const auto &b : this->bnds)
    {
        n += b.nSegments();
    }

    return n;
}

int SelectorMBond::nCutGroups() const
{
    int n = 0;

    for (const auto &b : this->bnds)
    {
        n += b.nCutGroups();
    }

    return n;
}

int SelectorMBond::nMolecules() const
{
    return this->bnds.count();
}

bool SelectorMBond::isEmpty() const
{
    return this->bnds.isEmpty();
}

SelectorMBond::const_iterator SelectorMBond::begin() const
{
    return this->bnds.constBegin();
}

SelectorMBond::const_iterator SelectorMBond::end() const
{
    return this->bnds.constEnd();
}

SelectorMBond::const_iterator SelectorMBond::constBegin() const
{
    return this->bnds.constBegin();
}

SelectorMBond::const_iterator SelectorMBond::constEnd() const
{
    return this->bnds.constEnd();
}

QString SelectorMBond::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("SelectorMBond::empty");
    }
    else
    {
        QStringList parts;

        const auto n = this->count();

        if (n <= 10)
        {
            for (int i=0; i<n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }
        }
        else
        {
            for (int i=0; i<5; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }

            parts.append("...");

            for (int i=n-5; i<n; ++i)
            {
                const auto view = this->operator[](i);

                parts.append(QString("%1: %2 %3")
                    .arg(i).arg(view.data().number().toString())
                    .arg(view.toString()));
            }
        }

        return QObject::tr("SelectorMBond( size=%2\n%3\n)")
                    .arg(n)
                    .arg(parts.join("\n"));
    }
}

