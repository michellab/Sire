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

#include "selectormimproper.h"

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

RegisterMetaType<SelectorMImproper> r_simp;

/** Serialise to a binary datastream */
SIREMOL_EXPORT QDataStream &operator<<(QDataStream &ds, const SelectorMImproper &imps)
{
    writeHeader(ds, r_simp, 1);

    SharedDataStream sds(ds);

    sds << imps.imps << static_cast<const Property&>(imps);

    return ds;
}

/** Extract from a binary datastream */
SIREMOL_EXPORT QDataStream &operator>>(QDataStream &ds, SelectorMImproper &imps)
{
    VersionID v = readHeader(ds, r_simp);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> imps.imps >> static_cast<Property&>(imps);
    }
    else
        throw version_error(v, "1", r_simp, CODELOC);

    return ds;
}

SelectorMImproper::SelectorMImproper()
                  : ConcreteProperty<SelectorMImproper, Property>()
{}

SelectorMImproper::SelectorMImproper(const Improper &view)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not view.isEmpty())
        imps.append(SelectorImproper(view));
}

SelectorMImproper::SelectorMImproper(const Molecules &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
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

        this->imps.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorImproper d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->imps.append(d);
        }
    }
}

SelectorMImproper::SelectorMImproper(const MoleculeGroup &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->imps.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorImproper d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->imps.append(d);
        }
    }
}

SelectorMImproper::SelectorMImproper(const MolGroupsBase &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not mols.isEmpty())
    {
        const auto molnums = mols.molNums();
        this->imps.reserve(molnums.count());

        for (const auto &molnum : molnums)
        {
            SelectorImproper d(mols.at(molnum), map);

            if (not d.isEmpty())
                this->imps.append(d);
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectResult &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not mols.isEmpty())
    {
        this->imps.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorImproper d;

            if (mol->isA<SelectorImproper>())
                d = mol->asA<SelectorImproper>();
            else
                d = SelectorImproper(*mol, map);

            if (not d.isEmpty())
                this->imps.append(d);
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectResult &mols,
                                     const ImproperID &improper,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not mols.isEmpty())
    {
        this->imps.reserve(mols.count());

        for (const auto &mol : mols)
        {
            try
            {
                auto d = SelectorImproper(*mol, improper, map);

                if (not d.isEmpty())
                    this->imps.append(d);
            }
            catch(...)
            {}
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectorImproper &impropers)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not impropers.isEmpty())
        imps.append(impropers);
}

SelectorMImproper::SelectorMImproper(const SelectorMol &mols,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    if (not mols.isEmpty())
    {
        this->imps.reserve(mols.count());

        for (const auto &mol : mols)
        {
            SelectorImproper d(mol, map);

            if (not d.isEmpty())
                imps.append(d);
        }
    }
}

void SelectorMImproper::_append(const SelectorImproper &impropers)
{
    if (impropers.isEmpty())
        return;

    if (this->imps.isEmpty())
        this->imps.append(impropers);
    else
    {
        for (int i=0; i<impropers.count(); ++i)
        {
            this->_append(impropers(i));
        }
    }
}

void SelectorMImproper::_append(const Improper &improper)
{
    if (this->imps.isEmpty())
    {
        this->imps.append(SelectorImproper(improper));
    }
    else if (this->imps.last().data().number() != improper.data().number())
    {
        // new molecule
        this->imps.append(SelectorImproper(improper));
    }
    else
    {
        // a new view in the current molecule
        this->imps.last() = this->imps.last().add(improper);
    }
}

SelectorMImproper::SelectorMImproper(const SelectorMImproper &impropers,
                                     const SireBase::Slice &slice)
                  : SireBase::ConcreteProperty<SelectorMImproper,Property>()
{
    for (auto it = slice.begin(impropers.count());
         not it.atEnd(); it.next())
    {
        this->_append(impropers[it.value()]);
    }
}

SelectorMImproper::SelectorMImproper(const SelectorMImproper &impropers,
                                     const QList<qint64> &idxs)
                  : SireBase::ConcreteProperty<SelectorMImproper,Property>()
{
    for (const auto &idx : idxs)
    {
        this->_append(impropers[idx]);
    }
}

SelectorMImproper::SelectorMImproper(const SelectorM<Atom> &atoms,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
{
    for (const auto &mol_atoms : atoms)
    {
        const auto impropers = SelectorImproper(mol_atoms, map);
        this->_append(impropers);
    }
}

SelectorMImproper::SelectorMImproper(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const PropertyMap &map)
              : ConcreteProperty<SelectorMImproper, Property>()
{
    for (const auto &mol_atoms0 : atoms0)
    {
        for (const auto &mol_atoms1 : atoms1)
        {
            if (mol_atoms0.isSameMolecule(mol_atoms1))
            {
                const auto impropers = SelectorImproper(mol_atoms0,
                                                        mol_atoms1, map);
                this->_append(impropers);
            }
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const SelectorM<Atom> &atoms2,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
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
                        const auto impropers = SelectorImproper(
                                                           mol_atoms0,
                                                           mol_atoms1,
                                                           mol_atoms2, map);
                        this->_append(impropers);
                    }
                }
            }
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectorM<Atom> &atoms0,
                                     const SelectorM<Atom> &atoms1,
                                     const SelectorM<Atom> &atoms2,
                                     const SelectorM<Atom> &atoms3,
                                     const PropertyMap &map)
                  : ConcreteProperty<SelectorMImproper, Property>()
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
                                const auto impropers = SelectorImproper(
                                                                   mol_atoms0,
                                                                   mol_atoms1,
                                                                   mol_atoms2,
                                                                   mol_atoms3,
                                                                   map);
                                this->_append(impropers);
                            }
                        }
                    }
                }
            }
        }
    }
}

SelectorMImproper::SelectorMImproper(const SelectorMImproper &other)
                  : ConcreteProperty<SelectorMImproper, Property>(),
                    imps(other.imps)
{}

SelectorMImproper::~SelectorMImproper()
{}

const char* SelectorMImproper::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorMImproper>());
}

SelectorMImproper& SelectorMImproper::operator=(const SelectorMImproper &other)
{
    if (this != &other)
    {
        imps = other.imps;
        Property::operator=(other);
    }

    return *this;
}

bool SelectorMImproper::operator==(const SelectorMImproper &other) const
{
    return imps == other.imps;
}

bool SelectorMImproper::operator!=(const SelectorMImproper &other) const
{
    return not operator==(other);
}

Improper SelectorMImproper::operator[](int i) const
{
    i = SireID::Index(i).map(this->count());

    for (const auto &d : imps)
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

    return Improper();
}

SelectorMImproper SelectorMImproper::operator[](const SireBase::Slice &slice) const
{
    return SelectorMImproper(*this, slice);
}

SelectorMImproper SelectorMImproper::operator[](const QList<qint64> &idxs) const
{
    return SelectorMImproper(*this, idxs);
}

SelectorMImproper SelectorMImproper::operator[](const ImproperID &id) const
{
    SelectorMImproper ret;

    for (const auto &d : imps)
    {
        try
        {
            auto r = d(id);

            if (not r.isEmpty())
            {
                ret.imps.append(r);
            }
        }
        catch(...)
        {}
    }

    return ret;
}

Improper SelectorMImproper::operator()(int i) const
{
    return this->operator[](i);
}

SelectorMImproper SelectorMImproper::operator()(const SireBase::Slice &slice) const
{
    return this->operator[](slice);
}

SelectorMImproper SelectorMImproper::operator()(const QList<qint64> &idxs) const
{
    return this->operator[](idxs);
}

SelectorMImproper SelectorMImproper::operator()(const ImproperID &id) const
{
    return this->operator[](id);
}

QList<MolViewPtr> SelectorMImproper::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(imps.count());

    for (const auto &imp : imps)
    {
        l.append(MolViewPtr(imp.clone()));
    }

    return l;
}

int SelectorMImproper::count() const
{
    int n = 0;

    for (const auto &d : imps)
    {
        n += d.count();
    }

    return n;
}

int SelectorMImproper::size() const
{
    return this->count();
}

EvaluatorM SelectorMImproper::evaluate() const
{
    return EvaluatorM(this->atoms());
}

MoleculeGroup SelectorMImproper::toMoleculeGroup() const
{
    MoleculeGroup grp;

    for (const auto &d : this->imps)
    {
        grp.add(d);
    }

    return grp;

}

SelectResult SelectorMImproper::toSelectResult() const
{
    QList<MolViewPtr> r;

    for (const auto &d : imps)
    {
        r.append(d);
    }

    return SelectResult(r);
}

Molecule SelectorMImproper::molecule(int i) const
{
    return this->molecules().molecule(i);
}

Molecule SelectorMImproper::molecule(const QString &name) const
{
    return this->molecules().molecule(name);
}

Molecule SelectorMImproper::molecule(const MolID &molid)
{
    return this->molecules().molecule(molid);
}

SelectorMol SelectorMImproper::molecules() const
{
    QList<Molecule> mols;

    for (const auto &d : this->imps)
    {
        mols.append(d.molecule());
    }

    return SelectorMol(mols);
}

SelectorMol SelectorMImproper::molecules(int i) const
{
    return this->molecules().molecules(i);
}

SelectorMol SelectorMImproper::molecules(const SireBase::Slice &slice) const
{
    return this->molecules().molecules(slice);
}

SelectorMol SelectorMImproper::molecules(const QList<qint64> &idxs) const
{
    return this->molecules().molecules(idxs);
}

SelectorMol SelectorMImproper::molecules(const QString &name) const
{
    return this->molecules().molecules(name);
}

SelectorMol SelectorMImproper::molecules(const MolID &molid) const
{
    return this->molecules().molecules(molid);
}

Atom SelectorMImproper::atom(int i) const
{
    return this->atoms()(i);
}

Atom SelectorMImproper::atom(const QString &name) const
{
    return this->atoms()(name);
}

Atom SelectorMImproper::atom(const AtomID &atomid) const
{
    return this->atoms()(atomid);
}

Residue SelectorMImproper::residue(int i) const
{
    return this->residues()(i);
}

Residue SelectorMImproper::residue(const QString &name) const
{
    return this->residues()(name);
}

Residue SelectorMImproper::residue(const ResID &resid) const
{
    return this->residues()(resid);
}

Chain SelectorMImproper::chain(int i) const
{
    return this->chains()(i);
}

Chain SelectorMImproper::chain(const QString &name) const
{
    return this->chains()(name);
}

Chain SelectorMImproper::chain(const ChainID &chainid) const
{
    return this->chains()(chainid);
}

Segment SelectorMImproper::segment(int i) const
{
    return this->segments()(i);
}

Segment SelectorMImproper::segment(const QString &name) const
{
    return this->segments()(name);
}

Segment SelectorMImproper::segment(const SegID &segid) const
{
    return this->segments()(segid);
}

CutGroup SelectorMImproper::cutGroup(int i) const
{
    return this->cutGroups()(i);
}

CutGroup SelectorMImproper::cutGroup(const QString &name) const
{
    return this->cutGroups()(name);
}

CutGroup SelectorMImproper::cutGroup(const CGID &cgid) const
{
    return this->cutGroups()(cgid);
}

SelectorM<Atom> SelectorMImproper::atoms() const
{
    QList< Selector<Atom> > ret;

    for (const auto &d : this->imps)
    {
        ret.append(d.atoms());
    }

    return SelectorM<Atom>(ret);
}

SelectorM<Atom> SelectorMImproper::atoms(int i) const
{
    return this->atoms().atoms(i);
}

SelectorM<Atom> SelectorMImproper::atoms(const SireBase::Slice &slice) const
{
    return this->atoms().atoms(slice);
}

SelectorM<Atom> SelectorMImproper::atoms(const QList<qint64> &idxs) const
{
    return this->atoms().atoms(idxs);
}

SelectorM<Atom> SelectorMImproper::atoms(const QString &name) const
{
    return this->atoms().atoms(name);
}

SelectorM<Atom> SelectorMImproper::atoms(const AtomID &atomid) const
{
    return this->atoms().atoms(atomid);
}

SelectorM<Residue> SelectorMImproper::residues() const
{
    QList< Selector<Residue> > ret;

    for (const auto &d : this->imps)
    {
        ret.append(d.residues());
    }

    return SelectorM<Residue>(ret);
}

SelectorM<Residue> SelectorMImproper::residues(int i) const
{
    return this->residues().residues(i);
}

SelectorM<Residue> SelectorMImproper::residues(const SireBase::Slice &slice) const
{
    return this->residues().residues(slice);
}

SelectorM<Residue> SelectorMImproper::residues(const QList<qint64> &idxs) const
{
    return this->residues().residues(idxs);
}

SelectorM<Residue> SelectorMImproper::residues(const QString &name) const
{
    return this->residues().residues(name);
}

SelectorM<Residue> SelectorMImproper::residues(const ResID &resid) const
{
    return this->residues().residues(resid);
}

SelectorM<Chain> SelectorMImproper::chains() const
{
    QList< Selector<Chain> > ret;

    for (const auto &d : this->imps)
    {
        ret.append(d.chains());
    }

    return SelectorM<Chain>(ret);
}

SelectorM<Chain> SelectorMImproper::chains(int i) const
{
    return this->chains().chains(i);
}

SelectorM<Chain> SelectorMImproper::chains(const SireBase::Slice &slice) const
{
    return this->chains().chains(slice);
}

SelectorM<Chain> SelectorMImproper::chains(const QList<qint64> &idxs) const
{
    return this->chains().chains(idxs);
}

SelectorM<Chain> SelectorMImproper::chains(const QString &name) const
{
    return this->chains().chains(name);
}

SelectorM<Chain> SelectorMImproper::chains(const ChainID &chainid) const
{
    return this->chains().chains(chainid);
}

SelectorM<Segment> SelectorMImproper::segments() const
{
    QList< Selector<Segment> > ret;

    for (const auto &d : this->imps)
    {
        ret.append(d.segments());
    }

    return SelectorM<Segment>(ret);
}

SelectorM<Segment> SelectorMImproper::segments(int i) const
{
    return this->segments().segments(i);
}

SelectorM<Segment> SelectorMImproper::segments(const SireBase::Slice &slice) const
{
    return this->segments().segments(slice);
}

SelectorM<Segment> SelectorMImproper::segments(const QList<qint64> &idxs) const
{
    return this->segments().segments(idxs);
}

SelectorM<Segment> SelectorMImproper::segments(const QString &name) const
{
    return this->segments().segments(name);
}

SelectorM<Segment> SelectorMImproper::segments(const SegID &segid) const
{
    return this->segments().segments(segid);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups() const
{
    QList< Selector<CutGroup> > ret;

    for (const auto &d : this->imps)
    {
        ret.append(d.cutGroups());
    }

    return SelectorM<CutGroup>(ret);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups(int i) const
{
    return this->cutGroups().cutGroups(i);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups(const SireBase::Slice &slice) const
{
    return this->cutGroups().cutGroups(slice);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups(const QList<qint64> &idxs) const
{
    return this->cutGroups().cutGroups(idxs);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups(const QString &name) const
{
    return this->cutGroups().cutGroups(name);
}

SelectorM<CutGroup> SelectorMImproper::cutGroups(const CGID &cgid) const
{
    return this->cutGroups().cutGroups(cgid);
}

SelectResult SelectorMImproper::search(const QString &search_string) const
{
    return this->toSelectResult().search(search_string);
}

QList<ImproperID> SelectorMImproper::IDs() const
{
    QList<ImproperID> ret;

    for (const auto &d : this->imps)
    {
        ret += d.IDs();
    }

    return ret;
}

int SelectorMImproper::nAtoms() const
{
    int n = 0;

    for (const auto &d : this->imps)
    {
        n += d.nAtoms();
    }

    return n;
}

int SelectorMImproper::nResidues() const
{
    int n = 0;

    for (const auto &d : this->imps)
    {
        n += d.nResidues();
    }

    return n;
}

int SelectorMImproper::nChains() const
{
    int n = 0;

    for (const auto &d : this->imps)
    {
        n += d.nChains();
    }

    return n;
}

int SelectorMImproper::nSegments() const
{
    int n = 0;

    for (const auto &d : this->imps)
    {
        n += d.nSegments();
    }

    return n;
}

int SelectorMImproper::nCutGroups() const
{
    int n = 0;

    for (const auto &d : this->imps)
    {
        n += d.nCutGroups();
    }

    return n;
}

int SelectorMImproper::nMolecules() const
{
    return this->imps.count();
}

bool SelectorMImproper::isEmpty() const
{
    return this->imps.isEmpty();
}

SelectorMImproper::const_iterator SelectorMImproper::begin() const
{
    return this->imps.constBegin();
}

SelectorMImproper::const_iterator SelectorMImproper::end() const
{
    return this->imps.constEnd();
}

SelectorMImproper::const_iterator SelectorMImproper::constBegin() const
{
    return this->imps.constBegin();
}

SelectorMImproper::const_iterator SelectorMImproper::constEnd() const
{
    return this->imps.constEnd();
}

QString SelectorMImproper::toString() const
{
    if (this->isEmpty())
    {
        return QObject::tr("SelectorMImproper::empty");
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

        return QObject::tr("SelectorMImproper( size=%2\n%3\n)")
                    .arg(n)
                    .arg(parts.join("\n"));
    }
}

SelectorMImproper SelectorMImproper::add(const SelectorMImproper &other) const
{
    SelectorMImproper ret(*this);

    for (const auto &value : other)
    {
        if (ret.isEmpty())
        {
            ret.imps.append(value);
        }
        else if (ret.imps.last().isSameMolecule(value))
        {
            for (int i=0; i<value.count(); ++i)
            {
                ret._append(value(i));
            }
        }
        else
        {
            ret.imps.append(value);
        }
    }

    return ret;
}

SelectorMImproper SelectorMImproper::intersection(const SelectorMImproper &other) const
{
    if (this->count() < other.count())
        return other.intersection(*this);

    SelectorMImproper ret;

    for (const auto &val : imps)
    {
        SelectorImproper intersect;

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

        ret.imps.append(intersect);
    }

    return ret;
}

SelectorMImproper SelectorMImproper::invert(const SireBase::PropertyMap &map) const
{
    SelectorMImproper ret;

    for (const auto &val : imps)
    {
        ret.imps.append(val.invert(map));
    }

    return ret;
}

SelectorMImproper SelectorMImproper::invert() const
{
    return this->invert(PropertyMap());
}

bool SelectorMImproper::hasProperty(const SireBase::PropertyName &key) const
{
    for (const auto &val : imps)
    {
        if (val.hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorMImproper::hasMetadata(const SireBase::PropertyName &key) const
{
    for (const auto &val : imps)
    {
        if (val.hasMetadata(key))
            return true;
    }

    return false;
}

bool SelectorMImproper::hasMetadata(const SireBase::PropertyName &key,
                                 const SireBase::PropertyName &metakey) const
{
    for (const auto &val : imps)
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

QStringList SelectorMImproper::propertyKeys() const
{
    QSet<QString> keys;

    for (const auto &val : imps)
    {
        keys += _to_set(val.propertyKeys());
    }

    return keys.values();
}

QStringList SelectorMImproper::metadataKeys() const
{
    QSet<QString> keys;

    for (const auto &val : imps)
    {
        keys += _to_set(val.metadataKeys());
    }

    return keys.values();
}

QStringList SelectorMImproper::metadataKeys(const SireBase::PropertyName &key) const
{
    QSet<QString> keys;

    for (const auto &val : imps)
    {
        keys += _to_set(val.metadataKeys(key));
    }

    return keys.values();
}

QList<SireBase::Properties> SelectorMImproper::properties() const
{
    QList<SireBase::Properties> props;

    for (const auto &val : imps)
    {
        props += val.properties();
    }

    return props;
}

QList<SireBase::PropertyPtr> SelectorMImproper::property(const SireBase::PropertyName &key) const
{
    QList<SireBase::PropertyPtr> props;

    bool has_prop = false;

    for (const auto &val : imps)
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
            "None of the impropers in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<SireBase::PropertyPtr> SelectorMImproper::property(const SireBase::PropertyName &key,
                                                      const Property &default_value) const
{
    QList<SireBase::PropertyPtr> props;

    for (const auto &val : imps)
    {
        props += val.property(key, default_value);
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorMImproper::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<SireUnits::Dimension::Angle> SelectorMImproper::sizes(const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> a;

    for (const auto &val : imps)
    {
        a += val.sizes(map);
    }

    return a;
}

QList<SireUnits::Dimension::Angle> SelectorMImproper::measures() const
{
    return this->sizes();
}

QList<SireUnits::Dimension::Angle> SelectorMImproper::measures(const SireBase::PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireCAS::Expression> SelectorMImproper::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<SireCAS::Expression> SelectorMImproper::potentials(const SireBase::PropertyMap &map) const
{
    QList<SireCAS::Expression> e;

    for (const auto &val : imps)
    {
        e += val.potentials(map);
    }

    return e;
}

QList<SireUnits::Dimension::MolarEnergy> SelectorMImproper::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::MolarEnergy> SelectorMImproper::energies(
                                    const SireBase::PropertyMap &map) const
{
    QList<SireUnits::Dimension::MolarEnergy> e;

    for (const auto &val : imps)
    {
        e += val.energies(map);
    }

    return e;
}

SireUnits::Dimension::MolarEnergy SelectorMImproper::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy SelectorMImproper::energy(
                                const SireBase::PropertyMap &map) const
{
    SireUnits::Dimension::MolarEnergy nrg(0);

    for (const auto &val : imps)
    {
        nrg += val.energy(map);
    }

    return nrg;
}
