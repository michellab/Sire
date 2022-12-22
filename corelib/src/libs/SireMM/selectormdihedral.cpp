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

#include "selectormdihedral.h"

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

RegisterMetaType<SelectorMDihedral> r_sdih;

/** Serialise to a binary datastream */
SIREMM_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMDihedral &dihs)
{
    writeHeader(ds, r_sdih, 1);

    SharedDataStream sds(ds);

    sds << dihs.dihs << static_cast<const Property&>(dihs);

    return ds;
}

/** Extract from a binary datastream */
SIREMM_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMDihedral &dihs)
{
    VersionID v = readHeader(ds, r_sdih);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> dihs.dihs >> static_cast<Property&>(dihs);
    }
    else
        throw version_error(v, "1", r_sdih, CODELOC);

    return ds;
}

SelectorMDihedral::SelectorMDihedral()
                  : ConcreteProperty<SelectorMDihedral, Property>()
{}

SelectorMDihedral::SelectorMDihedral(const Dihedral &view)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not view.isEmpty())
        dihs.append(SelectorDihedral(view));
}

SelectorMDihedral::SelectorMDihedral(const Molecules &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
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

        this->dihs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorDihedral d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->dihs.append(d);
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const MoleculeGroup &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->dihs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorDihedral d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->dihs.append(d);
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const MolGroupsBase &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->dihs.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorDihedral d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->dihs.append(d);
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectResult &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not mols.isEmpty())
    {
        this->dihs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorDihedral d;

            if (mol->isA<SelectorDihedral>())
                d = mol->asA<SelectorDihedral>();
            else
                d = SelectorDihedral(*mol, map);

            if (not d.isEmpty())
                this->dihs.append(d);
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectResult &mols,
                                     const DihedralID &dihedral,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not mols.isEmpty())
    {
        this->dihs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            try
            {
                auto d = SelectorDihedral(*mol, dihedral, map);

                if (not d.isEmpty())
                    this->dihs.append(d);
            }
            catch(...)
            {}
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorDihedral &dihedrals)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not dihedrals.isEmpty())
        dihs.append(dihedrals);
}

SelectorMDihedral::SelectorMDihedral(const SelectorMol &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    if (not mols.isEmpty())
    {
        this->dihs.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorDihedral d(mol, map);

            if (not d.isEmpty())
                dihs.append(d);
        }
    }
}

void SelectorMDihedral::_append(const SelectorDihedral &dihedrals)
{
    if (dihedrals.isEmpty())
        return;

    if (this->dihs.isEmpty())
        this->dihs.append(dihedrals);
    else
    {
        for (int i=0; i<dihedrals.count(); ++i)
        {
            this->_append(dihedrals(i));
        }
    }
}

void SelectorMDihedral::_append(const Dihedral &dihedral)
{
    if (this->dihs.isEmpty())
    {
        this->dihs.append(SelectorDihedral(dihedral));
    }
    else if (this->dihs.last().data().number() != dihedral.data().number())
    {
        // new molecule
        this->dihs.append(SelectorDihedral(dihedral));
    }
    else
    {
        // a new view in the current molecule
        this->dihs.last() = this->dihs.last().add(dihedral);
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorMDihedral &dihedrals,
                                     const SireBase::Slice &slice)
                  : SireBase::ConcreteProperty<SelectorMDihedral,Property>()
{
    for (auto it = slice.begin(dihedrals.count());
         not it.atEnd(); it.next())
    {
        this->_append(dihedrals[it.value()]);
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorMDihedral &dihedrals,
                                     const QList<qint64> &idxs)
                  : SireBase::ConcreteProperty<SelectorMDihedral,Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(dihedrals[idx]);
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorM<Atom> &atoms,
                               const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    for (const auto &mol_atoms : atoms)
    {
        const auto dihedrals = SelectorDihedral(mol_atoms, map);
        this->_append(dihedrals);
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const PropertyMap &map)
              : ConcreteProperty<SelectorMDihedral, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                const auto dihedrals = SelectorDihedral(mol_atoms0,
                                                        mol_atoms1, map);
                this->_append(dihedrals);
            }
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const SelectorM<Atom> &atoms2,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                for (const auto &mol_atoms2 : atoms2)
                {
                    if (mol_atoms0.isSameMolecule(mol_atoms2))
                    {
                        const auto dihedrals = SelectorDihedral(
                                                           mol_atoms0,
                                                           mol_atoms1,
                                                           mol_atoms2, map);
                        this->_append(dihedrals);
                    }
                }
            }
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const SelectorM<Atom> &atoms2,
                                     const SelectorM<Atom> &atoms3,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMDihedral, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                for (const auto &mol_atoms2 : atoms2)
                {
                    if (mol_atoms0.isSameMolecule(mol_atoms2))
                    {
                        for (const auto &mol_atoms3 : atoms3)
                        {
                            if (mol_atoms0.isSameMolecule(mol_atoms3))
                            {
                                const auto dihedrals = SelectorDihedral(
                                                                   mol_atoms0,
                                                                   mol_atoms1,
                                                                   mol_atoms2,
                                                                   mol_atoms3,
                                                                   map);
                                this->_append(dihedrals);
                            }
                        }
                    }
                }
            }
        }
    }
}

SelectorMDihedral::SelectorMDihedral(const SelectorMDihedral &other)
                  : ConcreteProperty<SelectorMDihedral, Property>(),
                    dihs(other.dihs)
{}

SelectorMDihedral::~SelectorMDihedral()
{}

const char* SelectorMDihedral::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMDihedral>());
}

SelectorMDihedral& SelectorMDihedral::operator=(const SelectorMDihedral &other)
{
    if (this != &other)
    {
        dihs = other.dihs;
        Property::operator=(other);
    }

    return *this;
}

bool SelectorMDihedral::operator==(const SelectorMDihedral &other) const
{
    return dihs == other.dihs;
}

bool SelectorMDihedral::operator!=(const SelectorMDihedral &other) const
{
    return not operator==(other);
}

Dihedral SelectorMDihedral::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    for (const auto &d : dihs)
    {
        if (i < d.count())
        {
            return d(i);
        }
        else
        {
            i -= d.count();
        }
    }

    throw SireError::program_bug(QObject::tr("Should not get here!"), CODELOC);

    return Dihedral();
}

SelectorMDihedral SelectorMDihedral::operator[](const SireBase::Slice &slice) const
{
    return SelectorMDihedral(*this, slice);
}

SelectorMDihedral SelectorMDihedral::operator[](const QList<qint64> &idxs) const
{
    return SelectorMDihedral(*this, idxs);
}

SelectorMDihedral SelectorMDihedral::operator[](const DihedralID &id) const
{
    SelectorMDihedral ret;

    for (const auto &d : dihs)
    {
        try
        {
            auto r = d(id);

            if (not r.isEmpty())
            {
                ret.dihs.append(r);
            }
        }
        catch(...)
        {}
    }

    return ret;
}

Dihedral SelectorMDihedral::operator()(int i) const
{
    return this->operator[](i);
}

SelectorMDihedral SelectorMDihedral::operator()(const SireBase::Slice &slice) const
{
    return this->operator[](slice);
}

SelectorMDihedral SelectorMDihedral::operator()(const QList<qint64> &idxs) const
{
    return this->operator[](idxs);
}

SelectorMDihedral SelectorMDihedral::operator()(const DihedralID &id) const
{
    return this->operator[](id);
}

bool SelectorMDihedral::isSelector() const
{
    return true;
}

QList<MolViewPtr> SelectorMDihedral::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(dihs.count());

    for (const auto &dih : dihs)
    {
        l.append(MolViewPtr(dih.clone()));
    }

    return l;
}

Molecules SelectorMDihedral::toMolecules() const
{
    return Molecules(this->dihs);
}

void SelectorMDihedral::update(const Molecules &molecules)
{
    // better to create a map from MolNum to index here
    QMultiHash<MolNum, int> molnum_to_idx;
    molnum_to_idx.reserve(this->dihs.count());

    int i = 0;

    for (const auto &mol : this->dihs)
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
            this->dihs[it.value()].update(mol.data());
        }
    }
}

int SelectorMDihedral::count() const
{
    int n = 0;

    for (const auto &d : dihs)
    {
        n += d.count();
    }

    return n;
}

int SelectorMDihedral::size() const
{
    return this->count();
}

EvaluatorM SelectorMDihedral::evaluate() const
{
    return EvaluatorM(this->atoms());
}

MoleculeGroup SelectorMDihedral::toMoleculeGroup() const
{
    MoleculeGroup grp;

    for (const auto &d : this->dihs)
    {
        grp.add(d);
    }

    return grp;

}

SelectResult SelectorMDihedral::toSelectResult() const
{
    QList<MolViewPtr> r;

    for (const auto &d : dihs)
    {
        r.append(d);
    }

    return SelectResult(r);
}

Molecule SelectorMDihedral::molecule(int i, const PropertyMap &map) const
{
    return this->molecules().molecule(i, map);
}

Molecule SelectorMDihedral::molecule(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecule(name, map);
}

Molecule SelectorMDihedral::molecule(const MolID &molid, const PropertyMap &map)
{
    return this->molecules().molecule(molid, map);
}

SelectorMol SelectorMDihedral::molecules() const
{
    QList<Molecule> mols;

    for (const auto &d : this->dihs)
    {
        mols.append(d.molecule());
    }

    return SelectorMol(mols);
}

SelectorMol SelectorMDihedral::molecules(int i, const PropertyMap &map) const
{
    return this->molecules().molecules(i, map);
}

SelectorMol SelectorMDihedral::molecules(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->molecules().molecules(slice, map);
}

SelectorMol SelectorMDihedral::molecules(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->molecules().molecules(idxs, map);
}

SelectorMol SelectorMDihedral::molecules(const QString &name, const PropertyMap &map) const
{
    return this->molecules().molecules(name, map);
}

SelectorMol SelectorMDihedral::molecules(const MolID &molid, const PropertyMap &map) const
{
    return this->molecules().molecules(molid, map);
}

Atom SelectorMDihedral::atom(int i, const PropertyMap &map) const
{
    return this->atoms().atom(i, map);
}

Atom SelectorMDihedral::atom(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atom(name, map);
}

Atom SelectorMDihedral::atom(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atom(atomid, map);
}

Residue SelectorMDihedral::residue(int i, const PropertyMap &map) const
{
    return this->residues().residue(i, map);
}

Residue SelectorMDihedral::residue(const QString &name, const PropertyMap &map) const
{
    return this->residues().residue(name, map);
}

Residue SelectorMDihedral::residue(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residue(resid, map);
}

Chain SelectorMDihedral::chain(int i, const PropertyMap &map) const
{
    return this->chains().chain(i, map);
}

Chain SelectorMDihedral::chain(const QString &name, const PropertyMap &map) const
{
    return this->chains().chain(name, map);
}

Chain SelectorMDihedral::chain(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chain(chainid, map);
}

Segment SelectorMDihedral::segment(int i, const PropertyMap &map) const
{
    return this->segments().segment(i, map);
}

Segment SelectorMDihedral::segment(const QString &name, const PropertyMap &map) const
{
    return this->segments().segment(name, map);
}

Segment SelectorMDihedral::segment(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segment(segid, map);
}

CutGroup SelectorMDihedral::cutGroup(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(i, map);
}

CutGroup SelectorMDihedral::cutGroup(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(name, map);
}

CutGroup SelectorMDihedral::cutGroup(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroup(cgid, map);
}

SelectorM<Atom> SelectorMDihedral::atoms() const
{
    QList< Selector<Atom> > ret;

    for (const auto &d : this->dihs)
    {
        ret.append(d.atoms());
    }

    return SelectorM<Atom>(ret);
}

SelectorM<Atom> SelectorMDihedral::atoms(int i, const PropertyMap &map) const
{
    return this->atoms().atoms(i, map);
}

SelectorM<Atom> SelectorMDihedral::atoms(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->atoms().atoms(slice, map);
}

SelectorM<Atom> SelectorMDihedral::atoms(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->atoms().atoms(idxs, map);
}

SelectorM<Atom> SelectorMDihedral::atoms(const QString &name, const PropertyMap &map) const
{
    return this->atoms().atoms(name, map);
}

SelectorM<Atom> SelectorMDihedral::atoms(const AtomID &atomid, const PropertyMap &map) const
{
    return this->atoms().atoms(atomid, map);
}

SelectorM<Residue> SelectorMDihedral::residues() const
{
    QList< Selector<Residue> > ret;

    for (const auto &d : this->dihs)
    {
        ret.append(d.residues());
    }

    return SelectorM<Residue>(ret);
}

SelectorM<Residue> SelectorMDihedral::residues(int i, const PropertyMap &map) const
{
    return this->residues().residues(i, map);
}

SelectorM<Residue> SelectorMDihedral::residues(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->residues().residues(slice, map);
}

SelectorM<Residue> SelectorMDihedral::residues(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->residues().residues(idxs, map);
}

SelectorM<Residue> SelectorMDihedral::residues(const QString &name, const PropertyMap &map) const
{
    return this->residues().residues(name, map);
}

SelectorM<Residue> SelectorMDihedral::residues(const ResID &resid, const PropertyMap &map) const
{
    return this->residues().residues(resid, map);
}

SelectorM<Chain> SelectorMDihedral::chains() const
{
    QList< Selector<Chain> > ret;

    for (const auto &d : this->dihs)
    {
        ret.append(d.chains());
    }

    return SelectorM<Chain>(ret);
}

SelectorM<Chain> SelectorMDihedral::chains(int i, const PropertyMap &map) const
{
    return this->chains().chains(i, map);
}

SelectorM<Chain> SelectorMDihedral::chains(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->chains().chains(slice, map);
}

SelectorM<Chain> SelectorMDihedral::chains(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->chains().chains(idxs, map);
}

SelectorM<Chain> SelectorMDihedral::chains(const QString &name, const PropertyMap &map) const
{
    return this->chains().chains(name, map);
}

SelectorM<Chain> SelectorMDihedral::chains(const ChainID &chainid, const PropertyMap &map) const
{
    return this->chains().chains(chainid, map);
}

SelectorM<Segment> SelectorMDihedral::segments() const
{
    QList< Selector<Segment> > ret;

    for (const auto &d : this->dihs)
    {
        ret.append(d.segments());
    }

    return SelectorM<Segment>(ret);
}

SelectorM<Segment> SelectorMDihedral::segments(int i, const PropertyMap &map) const
{
    return this->segments().segments(i, map);
}

SelectorM<Segment> SelectorMDihedral::segments(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->segments().segments(slice, map);
}

SelectorM<Segment> SelectorMDihedral::segments(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->segments().segments(idxs, map);
}

SelectorM<Segment> SelectorMDihedral::segments(const QString &name, const PropertyMap &map) const
{
    return this->segments().segments(name, map);
}

SelectorM<Segment> SelectorMDihedral::segments(const SegID &segid, const PropertyMap &map) const
{
    return this->segments().segments(segid, map);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups() const
{
    QList< Selector<CutGroup> > ret;

    for (const auto &d : this->dihs)
    {
        ret.append(d.cutGroups());
    }

    return SelectorM<CutGroup>(ret);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups(int i, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(i, map);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups(const SireBase::Slice &slice, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(slice, map);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups(const QList<qint64> &idxs, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(idxs, map);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups(const QString &name, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(name, map);
}

SelectorM<CutGroup> SelectorMDihedral::cutGroups(const CGID &cgid, const PropertyMap &map) const
{
    return this->cutGroups().cutGroups(cgid, map);
}


SelectResult SelectorMDihedral::search(const QString &search_string) const
{
    Select search(search_string);
    return search(this->toSelectResult());
}

SelectResult SelectorMDihedral::search(const QString &search_string,
                                       const PropertyMap &map) const
{
    Select search(search_string);
    return search(this->toSelectResult(), map);
}

QList<DihedralID> SelectorMDihedral::IDs() const
{
    QList<DihedralID> ret;

    for (const auto &d : this->dihs)
    {
        ret += d.IDs();
    }

    return ret;
}

int SelectorMDihedral::nAtoms() const
{
    int n = 0;

    for (const auto &d : this->dihs)
    {
        n += d.nAtoms();
    }

    return n;
}

int SelectorMDihedral::nResidues() const
{
    int n = 0;

    for (const auto &d : this->dihs)
    {
        n += d.nResidues();
    }

    return n;
}

int SelectorMDihedral::nChains() const
{
    int n = 0;

    for (const auto &d : this->dihs)
    {
        n += d.nChains();
    }

    return n;
}

int SelectorMDihedral::nSegments() const
{
    int n = 0;

    for (const auto &d : this->dihs)
    {
        n += d.nSegments();
    }

    return n;
}

int SelectorMDihedral::nCutGroups() const
{
    int n = 0;

    for (const auto &d : this->dihs)
    {
        n += d.nCutGroups();
    }

    return n;
}

int SelectorMDihedral::nMolecules() const
{
    return this->dihs.count();
}

bool SelectorMDihedral::isEmpty() const
{
    return this->dihs.isEmpty();
}

int SelectorMDihedral::nFrames() const
{
    return this->nFrames(PropertyMap());
}

int SelectorMDihedral::nFrames(const SireBase::PropertyMap &map) const
{
    return SireMol::detail::_nFrames(this->dihs, map);
}

void SelectorMDihedral::loadFrame(int frame)
{
    this->loadFrame(frame, PropertyMap());
}

void SelectorMDihedral::saveFrame(int frame)
{
    this->saveFrame(frame, PropertyMap());
}

void SelectorMDihedral::saveFrame()
{
    this->saveFrame(PropertyMap());
}

void SelectorMDihedral::deleteFrame(int frame)
{
    this->deleteFrame(frame, PropertyMap());
}

void SelectorMDihedral::loadFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_loadFrame(this->dihs, frame, map);
}

void SelectorMDihedral::saveFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->dihs, frame, map);
}

void SelectorMDihedral::saveFrame(const SireBase::PropertyMap &map)
{
    SireMol::detail::_saveFrame(this->dihs, map);
}

void SelectorMDihedral::deleteFrame(int frame, const SireBase::PropertyMap &map)
{
    SireMol::detail::_deleteFrame(this->dihs, frame, map);
}

SelectorMDihedral::const_iterator SelectorMDihedral::begin() const
{
    return this->dihs.constBegin();
}

SelectorMDihedral::const_iterator SelectorMDihedral::end() const
{
    return this->dihs.constEnd();
}

SelectorMDihedral::const_iterator SelectorMDihedral::constBegin() const
{
    return this->dihs.constBegin();
}

SelectorMDihedral::const_iterator SelectorMDihedral::constEnd() const
{
    return this->dihs.constEnd();
}

QString SelectorMDihedral::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("SelectorMDihedral::empty");
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

        return QObject::tr("SelectorMDihedral( size=%2\n%3\n)")
                    .arg(n)
                    .arg(parts.join("\n"));
    }
}

SelectorMDihedral SelectorMDihedral::add(const SelectorMDihedral &other) const
{
    SelectorMDihedral ret(*this);

    for (const auto &value : other)
    {
        if (ret.isEmpty())
        {
            ret.dihs.append(value);
        }
        else if (ret.dihs.last().isSameMolecule(value))
        {
            for (int i=0; i<value.count(); ++i)
            {
                ret._append(value(i));
            }
        }
        else
        {
            ret.dihs.append(value);
        }
    }

    return ret;
}

SelectorMDihedral SelectorMDihedral::intersection(const SelectorMDihedral &other) const
{
    if (this->count() < other.count())
        return other.intersection(*this);

    SelectorMDihedral ret;

    for (const auto &val : dihs)
    {
        SelectorDihedral intersect;

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

        ret.dihs.append(intersect);
    }

    return ret;
}

SelectorMDihedral SelectorMDihedral::invert(const SireBase::PropertyMap &map) const
{
    SelectorMDihedral ret;

    for (const auto &val : dihs)
    {
        ret.dihs.append(val.invert(map));
    }

    return ret;
}

SelectorMDihedral SelectorMDihedral::invert() const
{
    return this->invert(PropertyMap());
}

bool SelectorMDihedral::hasProperty(const SireBase::PropertyName &key) const
{
    for (const auto &val : dihs)
    {
        if (val.hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorMDihedral::hasMetadata(const SireBase::PropertyName &key) const
{
    for (const auto &val : dihs)
    {
        if (val.hasMetadata(key))
            return true;
    }

    return false;
}

bool SelectorMDihedral::hasMetadata(const SireBase::PropertyName &key,
                                 const SireBase::PropertyName &metakey) const
{
    for (const auto &val : dihs)
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

QStringList SelectorMDihedral::propertyKeys() const
{
    QSet<QString> keys;

    for (const auto &val : dihs)
    {
        keys += _to_set(val.propertyKeys());
    }

    return keys.values();
}

QStringList SelectorMDihedral::metadataKeys() const
{
    QSet<QString> keys;

    for (const auto &val : dihs)
    {
        keys += _to_set(val.metadataKeys());
    }

    return keys.values();
}

QStringList SelectorMDihedral::metadataKeys(const SireBase::PropertyName &key) const
{
    QSet<QString> keys;

    for (const auto &val : dihs)
    {
        keys += _to_set(val.metadataKeys(key));
    }

    return keys.values();
}

QList<SireBase::Properties> SelectorMDihedral::properties() const
{
    QList<SireBase::Properties> props;

    for (const auto &val : dihs)
    {
        props += val.properties();
    }

    return props;
}

QList<SireBase::PropertyPtr> SelectorMDihedral::property(const SireBase::PropertyName &key) const
{
    QList<SireBase::PropertyPtr> props;

    bool has_prop = false;

    for (const auto &val : dihs)
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
            "None of the dihedrals in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<SireBase::PropertyPtr> SelectorMDihedral::property(const SireBase::PropertyName &key,
                                                      const Property &default_value) const
{
    QList<SireBase::PropertyPtr> props;

    for (const auto &val : dihs)
    {
        props += val.property(key, default_value);
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorMDihedral::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<SireUnits::Dimension::Angle> SelectorMDihedral::sizes(const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> a;

    for (const auto &val : dihs)
    {
        a += val.sizes(map);
    }

    return a;
}

QList<SireUnits::Dimension::Angle> SelectorMDihedral::measures() const
{
    return this->sizes();
}

QList<SireUnits::Dimension::Angle> SelectorMDihedral::measures(const SireBase::PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireCAS::Expression> SelectorMDihedral::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<SireCAS::Expression> SelectorMDihedral::potentials(const SireBase::PropertyMap &map) const
{
    QList<SireCAS::Expression> e;

    for (const auto &val : dihs)
    {
        e += val.potentials(map);
    }

    return e;
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMDihedral::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::GeneralUnit> SelectorMDihedral::energies(
                                    const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::GeneralUnit> e;

    for (const auto &val : dihs)
    {
        e += val.energies(map);
    }

    return e;
}

SireUnits::Dimension::GeneralUnit SelectorMDihedral::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::GeneralUnit SelectorMDihedral::energy(
                                const SireBase::PropertyMap &map) const
{
    SireUnits::Dimension::GeneralUnit nrg(0);

    for (const auto &val : dihs)
    {
        nrg += val.energy(map);
    }

    return nrg;
}
