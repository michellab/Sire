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

#include "SireCAS/expression.h"

#include "SireMol/errors.h"
#include "SireError/errors.h"
#include "SireBase/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;
using namespace SireMM;
using namespace SireID;

RegisterMetaType<SelectorMBond> r_sbnd;

/** Serialise to a binary datastream */
SIREMM_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMBond &bnds)
{
    writeHeader(ds, r_sbnd, 1);

    SharedDataStream sds(ds);

    sds << bnds.bnds << static_cast<const Property&>(bnds);

    return ds;
}

/** Extract from a binary datastream */
SIREMM_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMBond &bnds)
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
    if (not view.isEmpty())
        bnds.append(SelectorBond(view));
}

SelectorMBond::SelectorMBond(const Molecules &mols,
                             const PropertyMap &map)
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
            SelectorBond b(mols.at(molnum), map);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MoleculeGroup &mols,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(mols.at(molnum), map);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const MolGroupsBase &mols,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->bnds.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorBond b(mols.at(molnum), map);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectResult &mols,
                             const PropertyMap &map)
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
                b = SelectorBond(*mol, map);

            if (not b.isEmpty())
                this->bnds.append(b);
        }
    }
}

SelectorMBond::SelectorMBond(const SelectResult &mols, const BondID &bond,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            try
            {
                auto b = SelectorBond(*mol, bond, map);

                if (not b.isEmpty())
                    this->bnds.append(b);
            }
            catch(...)
            {}
        }
    }
}

SelectorMBond::SelectorMBond(const SelectorM<Atom> &atoms,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    for (const auto &mol_atoms : atoms)
    {
        const auto bonds = SelectorBond(mol_atoms, map);
        this->_append(bonds);
    }
}

SelectorMBond::SelectorMBond(const SelectorM<Atom> &atoms0,
                             const SelectorM<Atom> &atoms1,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                const auto bonds = SelectorBond(mol_atoms0, mol_atoms1, map);
                this->_append(bonds);
            }
        }
    }
}

SelectorMBond::SelectorMBond(const SelectorBond &bonds)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not bonds.isEmpty())
        bnds.append(bonds);
}

SelectorMBond::SelectorMBond(const SelectorMol &mols,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorMBond, Property>()
{
    if (not mols.isEmpty())
    {
        this->bnds.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorBond b(mol, map);

            if (not b.isEmpty())
                bnds.append(b);
        }
    }
}

void SelectorMBond::_append(const SelectorBond &bonds)
{
    if (bonds.isEmpty())
        return;

    if (this->bnds.isEmpty())
        this->bnds.append(bonds);
    else
    {
        for (int i=0; i<bonds.count(); ++i)
        {
            this->_append(bonds(i));
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

QList<MolViewPtr> SelectorMBond::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(bnds.count());

    for (const auto &bnd : bnds)
    {
        l.append(MolViewPtr(bnd.clone()));
    }

    return l;
}

Molecules SelectorMBond::toMolecules() const
{
    return Molecules(this->bnds);
}

void SelectorMBond::update(const Molecules &molecules)
{
    // better to create a map from MolNum to index here
    QMultiHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(this->bnds.count());

    int i = 0;

    for (const auto &mol : this->bnds)
    {
        molnum_to_idx.insert(mol.data().number(), i);
        i += 1;
    }

    for (const auto &mol : molecules)
    {
        const auto molnum = mol.data().number();

        auto it = molnum_to_idx.constFind(molnum);

        while (it != molnum_to_idx.constEnd() && it.key() == molnum)
        {
            this->bnds[it.value()].update(mol.data());
        }
    }
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

bool SelectorMBond::isSelector() const
{
    return true;
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

Molecule SelectorMBond::molecule(int i, const PropertyMap &map) const
{
    return this->molecules().molecule(i, map);
}

Molecule SelectorMBond::molecule(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecule(name, map);
}

Molecule SelectorMBond::molecule(const MolID &molid, const PropertyMap &map)
{
    return this->molecules().molecule(molid, map);
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

SelectorMol SelectorMBond::molecules(int i, const PropertyMap &map) const
{
    return this->molecules().molecules(i, map);
}

SelectorMol SelectorMBond::molecules(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->molecules().molecules(slice, map);
}

SelectorMol SelectorMBond::molecules(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->molecules().molecules(idxs, map);
}

SelectorMol SelectorMBond::molecules(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecules(name, map);
}

SelectorMol SelectorMBond::molecules(const MolID &molid, const PropertyMap &map) const
{
    return this->molecules().molecules(molid, map);
}

Atom SelectorMBond::atom(int i, const PropertyMap &map) const
{
    return this->atoms().atom(i, map);
}

Atom SelectorMBond::atom(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atom(name, map);
}

Atom SelectorMBond::atom(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atom(atomid, map);
}

Residue SelectorMBond::residue(int i, const PropertyMap &map) const
{
    return this->residues().residue(i, map);
}

Residue SelectorMBond::residue(const QString &name, const PropertyMap &map) const
{
    return this->residues().residue(name, map);
}

Residue SelectorMBond::residue(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residue(resid, map);
}

Chain SelectorMBond::chain(int i, const PropertyMap &map) const
{
    return this->chains().chain(i, map);
}

Chain SelectorMBond::chain(const QString &name, const PropertyMap &map) const
{
    return this->chains().chain(name, map);
}

Chain SelectorMBond::chain(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chain(chainid, map);
}

Segment SelectorMBond::segment(int i, const PropertyMap &map) const
{
    return this->segments().segment(i, map);
}

Segment SelectorMBond::segment(const QString &name, const PropertyMap &map) const
{
    return this->segments().segment(name, map);
}

Segment SelectorMBond::segment(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segment(segid, map);
}

CutGroup SelectorMBond::cutGroup(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(i, map);
}

CutGroup SelectorMBond::cutGroup(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(name, map);
}

CutGroup SelectorMBond::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(cgid, map);
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

SelectorM<Atom> SelectorMBond::atoms(int i, const PropertyMap &map) const
{
    return this->atoms().atoms(i, map);
}

SelectorM<Atom> SelectorMBond::atoms(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->atoms().atoms(slice, map);
}

SelectorM<Atom> SelectorMBond::atoms(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->atoms().atoms(idxs, map);
}

SelectorM<Atom> SelectorMBond::atoms(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atoms(name, map);
}

SelectorM<Atom> SelectorMBond::atoms(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atoms(atomid, map);
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

SelectorM<Residue> SelectorMBond::residues(int i, const PropertyMap &map) const
{
    return this->residues().residues(i, map);
}

SelectorM<Residue> SelectorMBond::residues(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->residues().residues(slice, map);
}

SelectorM<Residue> SelectorMBond::residues(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->residues().residues(idxs, map);
}

SelectorM<Residue> SelectorMBond::residues(const QString &name, const PropertyMap &map) const
{
    return this->residues().residues(name, map);
}

SelectorM<Residue> SelectorMBond::residues(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residues(resid, map);
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

SelectorM<Chain> SelectorMBond::chains(int i, const PropertyMap &map) const
{
    return this->chains().chains(i, map);
}

SelectorM<Chain> SelectorMBond::chains(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->chains().chains(slice, map);
}

SelectorM<Chain> SelectorMBond::chains(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->chains().chains(idxs, map);
}

SelectorM<Chain> SelectorMBond::chains(const QString &name, const PropertyMap &map) const
{
    return this->chains().chains(name, map);
}

SelectorM<Chain> SelectorMBond::chains(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chains(chainid, map);
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

SelectorM<Segment> SelectorMBond::segments(int i, const PropertyMap &map) const
{
    return this->segments().segments(i, map);
}

SelectorM<Segment> SelectorMBond::segments(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->segments().segments(slice, map);
}

SelectorM<Segment> SelectorMBond::segments(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->segments().segments(idxs, map);
}

SelectorM<Segment> SelectorMBond::segments(const QString &name, const PropertyMap &map) const
{
    return this->segments().segments(name, map);
}

SelectorM<Segment> SelectorMBond::segments(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segments(segid, map);
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

SelectorM<CutGroup> SelectorMBond::cutGroups(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(i, map);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(slice, map);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(idxs, map);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(name, map);
}

SelectorM<CutGroup> SelectorMBond::cutGroups(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(cgid, map);
}

SelectResult SelectorMBond::search(const QString &search_string) const
{
    Select search(search_string);
    return search(this->toSelectResult());
}

SelectResult SelectorMBond::search(const QString &search_string,
                                   const PropertyMap &map) const
{
    Select search(search_string);
    return search(this->toSelectResult(), map);
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

int SelectorMBond::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int SelectorMBond::nFrames(const SireBase::PropertyMap &map) const
{
    return SireMol::detail::_nFrames(this->bnds, map);
}

void SelectorMBond::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void SelectorMBond::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void SelectorMBond::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void SelectorMBond::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void SelectorMBond::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_loadFrame(this->bnds, frame, map);
}

void SelectorMBond::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->bnds, frame, map);
}

void SelectorMBond::saveFrame(const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->bnds, map);
}

void SelectorMBond::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteFrame(this->bnds, frame, map);
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

SelectorMBond SelectorMBond::add(const SelectorMBond &other) const
{
    SelectorMBond ret(*this);

    for (const auto &value : other)
    {
        if (ret.isEmpty())
        {
            ret.bnds.append(value);
        }
        else if (ret.bnds.last().isSameMolecule(value))
        {
            for (int i=0; i<value.count(); ++i)
            {
                ret._append(value(i));
            }
        }
        else
        {
            ret.bnds.append(value);
        }
    }

    return ret;
}

SelectorMBond SelectorMBond::intersection(const SelectorMBond &other) const
{
    if (this->count() < other.count())
        return other.intersection(*this);

    SelectorMBond ret;

    for (const auto &val : bnds)
    {
        SelectorBond intersect;

        for (const auto &other_val : other)
        {
            if (val.isSameMolecule(other_val))
            {
                if (intersect.isEmpty())
                    intersect = val.intersection(other_val);
                else
                    intersect = intersect.add(val.intersection(other_val));
            }
        }

        ret.bnds.append(intersect);
    }

    return ret;
}

SelectorMBond SelectorMBond::invert(const SireBase::PropertyMap &map) const
{
    SelectorMBond ret;

    for (const auto &val : bnds)
    {
        ret.bnds.append(val.invert(map));
    }

    return ret;
}

SelectorMBond SelectorMBond::invert() const
{
    return this->invert(PropertyMap());
}

bool SelectorMBond::hasProperty(const SireBase::PropertyName &key) const
{
    for (const auto &val : bnds)
    {
        if (val.hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorMBond::hasMetadata(const SireBase::PropertyName &key) const
{
    for (const auto &val : bnds)
    {
        if (val.hasMetadata(key))
            return true;
    }

    return false;
}

bool SelectorMBond::hasMetadata(const SireBase::PropertyName &key,
                                const SireBase::PropertyName &metakey) const
{
    for (const auto &val : bnds)
    {
        if (val.hasMetadata(key, metakey))
            return true;
    }

    return false;
}

template<class T>
inline QSet<T> _to_set(const QList<T> &l)
{
    return QSet<T>(l.constBegin(), l.constEnd());
}

QStringList SelectorMBond::propertyKeys() const
{
    QSet<QString> keys;

    for (const auto &val : bnds)
    {
        keys += _to_set(val.propertyKeys());
    }

    return keys.values();
}

QStringList SelectorMBond::metadataKeys() const
{
    QSet<QString> keys;

    for (const auto &val : bnds)
    {
        keys += _to_set(val.metadataKeys());
    }

    return keys.values();
}

QStringList SelectorMBond::metadataKeys(const SireBase::PropertyName &key) const
{
    QSet<QString> keys;

    for (const auto &val : bnds)
    {
        keys += _to_set(val.metadataKeys(key));
    }

    return keys.values();
}

QList<SireBase::Properties> SelectorMBond::properties() const
{
    QList<SireBase::Properties> props;

    for (const auto &val : bnds)
    {
        props += val.properties();
    }

    return props;
}

QList<SireBase::PropertyPtr> SelectorMBond::property(const SireBase::PropertyName &key) const
{
    QList<SireBase::PropertyPtr> props;

    bool has_prop = false;

    for (const auto &val : bnds)
    {
        try
        {
            props += val.property(key);
            has_prop = true;
        }
        catch(const SireError::exception&)
        {
            PropertyPtr null(new NullProperty());

            for (int i=0; i<val.count(); ++i)
            {
                props.append(null);
            }
        }
    }

    if (not has_prop)
        throw SireBase::missing_property(QObject::tr(
            "None of the bonds in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<SireBase::PropertyPtr> SelectorMBond::property(const SireBase::PropertyName &key,
                                                     const Property &default_value) const
{
    QList<SireBase::PropertyPtr> props;

    for (const auto &val : bnds)
    {
        props += val.property(key, default_value);
    }

    return props;
}

QList<SireUnits::Dimension::Length> SelectorMBond::lengths() const
{
    return this->lengths(PropertyMap());
}

QList<SireUnits::Dimension::Length> SelectorMBond::lengths(const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::Length> l;

    for (const auto &val : bnds)
    {
        l += val.lengths(map);
    }

    return l;
}

QList<SireUnits::Dimension::Length> SelectorMBond::measures() const
{
    return this->lengths();
}

QList<SireUnits::Dimension::Length> SelectorMBond::measures(const SireBase::PropertyMap &map) const
{
    return this->lengths(map);
}

QList<SireCAS::Expression> SelectorMBond::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<SireCAS::Expression> SelectorMBond::potentials(const SireBase::PropertyMap &map) const
{
    QList<SireCAS::Expression> e;

    for (const auto &val : bnds)
    {
        e += val.potentials(map);
    }

    return e;
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMBond::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMBond::energies(
                                    const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::GeneralUnit> e;

    for (const auto &val : bnds)
    {
        e += val.energies(map);
    }

    return e;
}

SireUnits::Dimension::GeneralUnit SelectorMBond::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::GeneralUnit SelectorMBond::energy(
                                const SireBase::PropertyMap &map) const
{
    SireUnits::Dimension::GeneralUnit nrg(0);

    for (const auto &val : bnds)
    {
        nrg += val.energy(map);
    }

    return nrg;
}
