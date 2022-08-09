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

#include "selectorbond.h"

#include "twoatomfunctions.h"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"

#include "SireCAS/symbol.h"
#include "SireCAS/values.h"

#include "SireID/index.h"

#include "SireBase/slice.h"
#include "SireBase/errors.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireID;
using namespace SireStream;
using namespace SireError;
using namespace SireUnits;
using namespace SireUnits::Dimension;

using SireMM::detail::IDPair;

static const RegisterMetaType<SelectorBond> r_sbond;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const SelectorBond &sbond)
{
    writeHeader(ds, r_sbond, 1);
    SharedDataStream sds(ds);
    sds << sbond.bnds << static_cast<const MoleculeView&>(sbond);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, SelectorBond &sbond)
{
    auto v = readHeader(ds, r_sbond);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> sbond.bnds >> static_cast<MoleculeView&>(sbond);
    }
    else
        throw SireStream::version_error(v, "1", r_sbond, CODELOC);

    return ds;
}

SelectorBond::SelectorBond() : ConcreteProperty<SelectorBond, MoleculeView>()
{}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto bonds = c.getBonds();

        if (mol.selectedAll())
        {
            bnds = bonds;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &bond : bonds)
            {
                if (s.selected(bond.atom0()) and s.selected(bond.atom1()))
                {
                    bnds.append(bond);
                }
            }
        }
    }
}

SelectorBond::SelectorBond(const MoleculeData &moldata,
                           const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorBond, MoleculeView>()
{
    this->operator=(SelectorBond(Molecule(moldata), map));
}

inline
QSet<IDPair> _to_int_set(const QList<BondID> &vals,
                         const MoleculeInfoData &molinfo)
{
    QSet<IDPair> s;
    s.reserve(vals.count());

    for (const auto &val : vals)
    {
        s.insert(IDPair(molinfo.atomIdx(val[0]).value(),
                        molinfo.atomIdx(val[1]).value()));
    }

    return s;
}

inline
QList<BondID> _from_int_set(const QSet<IDPair> &vals)
{
    QVector<IDPair> v;
    v.reserve(vals.count());

    for (const auto &val : vals)
    {
        v.append(val);
    }

    std::sort(v.begin(), v.end());

    QList<BondID> l;
    l.reserve(v.count());

    for (const auto &val : v)
    {
        l.append(BondID(AtomIdx(val.atom0),
                        AtomIdx(val.atom1)));
    }

    return l;
}

inline
QList<AtomIdx> _filter(const QList<AtomIdx> &atoms, const AtomSelection &selection)
{
    QList<AtomIdx> ret;
    ret.reserve(atoms.count());

    for (const auto &atom : atoms)
    {
        if (selection.selected(atom))
        {
            ret.append(atom);
        }
    }

    return ret;
}

inline
QSet<IDPair> _filter(const QSet<IDPair> &bonds,
                     const QVector<quint32> &atoms,
                     int position)
{
    QSet<IDPair> result;

    result.reserve(bonds.count());

    for (const auto &bond : bonds)
    {
        const auto atomidx = bond[position];

        for (const auto &atom : atoms)
        {
            if (atomidx == atom)
            {
                result.insert(bond);
                break;
            }
        }
    }

    return result;
}

inline
QVector<quint32> _to_int(const Selector<Atom> &atoms)
{
    QVector<quint32> ret;

    const int n = atoms.count();

    ret.reserve(n);

    for (int i=0; i<n; ++i)
    {
        ret.append(atoms.index(i));
    }

    return ret;
}

inline
QVector<quint32> _to_int(const QList<AtomIdx> &atoms)
{
    QVector<quint32> ret;
    ret.reserve(atoms.count());

    for (const auto &atom : atoms)
    {
        ret.append(atom.value());
    }

    return ret;
}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const BondID &bond, const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(bond.atom0());
    auto atoms1 = mol.data().info().map(bond.atom1());

    if (mol.data().hasProperty(map["connectivity"]))
    {
        if (not mol.selectedAll())
        {
            const auto selection = mol.selection();
            atoms0 = _filter(atoms0, selection);
            atoms1 = _filter(atoms1, selection);
        }

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);

        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto bonds = _to_int_set(c.getBonds(), this->data().info());

        auto bonds01 = _filter(bonds, int_atoms0, 0);
        bonds01 = _filter(bonds01, int_atoms1, 1);

        auto bonds10 = _filter(bonds, int_atoms0, 1);
        bonds10 = _filter(bonds10, int_atoms1, 0);

        bnds = _from_int_set(bonds01 + bonds10);
    }
}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const QList<BondID> &bonds,
                           const SireBase::PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    if (bonds.count() == 1)
    {
        this->operator=(SelectorBond(mol, bonds[0], map));
        return;
    }
    else if (mol.data().hasProperty(map["connectivity"]) and not bonds.isEmpty())
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto all_bonds = _to_int_set(c.getBonds(), this->data().info());

        QSet<IDPair> selected_bonds;
        selected_bonds.reserve(all_bonds.count());

        for (const auto &bond : bonds)
        {
            auto atoms0 = mol.data().info().map(bond.atom0());
            auto atoms1 = mol.data().info().map(bond.atom1());

            if (not mol.selectedAll())
            {
                const auto selection = mol.selection();
                atoms0 = _filter(atoms0, selection);
                atoms1 = _filter(atoms1, selection);
            }

            auto int_atoms0 = _to_int(atoms0);
            auto int_atoms1 = _to_int(atoms1);

            auto bonds01 = _filter(all_bonds, int_atoms0, 0);
            bonds01 = _filter(bonds01, int_atoms1, 1);

            auto bonds10 = _filter(all_bonds, int_atoms0, 1);
            bonds10 = _filter(bonds10, int_atoms1, 0);

            selected_bonds += bonds01;
            selected_bonds += bonds10;

            if (selected_bonds.count() == all_bonds.count())
                break;
        }

        bnds = _from_int_set(selected_bonds);
    }
}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const AtomID &atom, const PropertyMap &map)
              : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    auto atoms = mol.data().info().map(atom);

    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        if (not mol.selectedAll())
        {
            const auto selection = mol.selection();
            atoms = _filter(atoms, selection);
        }

        auto int_atoms = _to_int(atoms);

        auto bonds = _to_int_set(c.getBonds(), this->data().info());

        auto bonds0 = _filter(bonds, int_atoms, 0);
        auto bonds1 = _filter(bonds, int_atoms, 1);
        auto bonds2 = _filter(bonds, int_atoms, 2);

        bnds = _from_int_set(bonds0 + bonds1 + bonds2);
    }
}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const AtomID &atom0, const AtomID &atom1,
                           const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    this->operator=(SelectorBond(mol, BondID(atom0, atom1), map));
}

SelectorBond::SelectorBond(const MoleculeData &mol,
                           const AtomID &atom, const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>()
{
    this->operator=(SelectorBond(Molecule(mol), atom, map));
}

SelectorBond::SelectorBond(const MoleculeData &mol,
                           const AtomID &atom0, const AtomID &atom1,
                           const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>()
{
    this->operator=(SelectorBond(Molecule(mol), atom0, atom1, map));
}

SelectorBond::SelectorBond(const Selector<Atom> &atoms,
                           const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>(atoms)
{
    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms = _to_int(atoms);

        auto bonds = _to_int_set(c.getBonds(), this->data().info());

        auto bonds0 = _filter(bonds, int_atoms, 0);
        auto bonds1 = _filter(bonds, int_atoms, 1);

        bnds = _from_int_set(bonds0 + bonds1);
    }
}

SelectorBond::SelectorBond(const Selector<Atom> &atoms0,
                           const Selector<Atom> &atoms1,
                           const PropertyMap &map)
             : ConcreteProperty<SelectorBond, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a bond from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);

        auto bonds = _to_int_set(c.getBonds(), this->data().info());

        auto bonds01 = _filter(bonds, int_atoms0, 0);
        bonds01 = _filter(bonds01, int_atoms1, 1);

        auto bonds10 = _filter(bonds, int_atoms0, 1);
        bonds10 = _filter(bonds10, int_atoms1, 0);

        bnds = _from_int_set(bonds01 + bonds10);
    }
}

SelectorBond::SelectorBond(const SelectorBond &other)
     : ConcreteProperty<SelectorBond, MoleculeView>(other), bnds(other.bnds)
{}

SelectorBond::~SelectorBond()
{}

const char* SelectorBond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorBond>());
}

SelectorBond& SelectorBond::operator=(const SelectorBond &other)
{
    if (this != &other)
    {
        bnds = other.bnds;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorBond::operator==(const SelectorBond &other) const
{
    return bnds == other.bnds and MoleculeView::operator==(other);
}

bool SelectorBond::operator!=(const SelectorBond &other) const
{
    return not operator==(other);
}

int SelectorBond::count() const
{
    return this->bnds.count();
}

int SelectorBond::size() const
{
    return this->count();
}

int SelectorBond::nViews() const
{
    return this->count();
}

QString SelectorBond::toString() const
{
    if (this->isNull())
        return QObject::tr("SelectorBond::null");
    else if (this->isEmpty())
        return QObject::tr("SelectorBond::empty");

    QStringList parts;

    if (this->count() <= 10)
    {
        for (int i=0; i<this->count(); ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i)
                            .arg(this->operator()(i).toString()));
        }
    }
    else
    {
        for (int i=0; i<5; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i)
                            .arg(this->operator()(i).toString()));
        }

        parts.append("...");

        for (int i=this->count()-5; i<this->count(); ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i)
                            .arg(this->operator()(i).toString()));
        }
    }

    return QObject::tr("SelectorBond( size=%1\n%2\n)")
                .arg(this->count()).arg(parts.join("\n"));
}

SelectorBond SelectorBond::add(const Bond &bond) const
{
    if (bond.isNull())
        return *this;

    if (this->isEmpty())
    {
        return SelectorBond(bond);
    }

    if (bond.data().number() != this->data().number())
    {
        throw SireError::incompatible_error(QObject::tr(
            "You cannot add bonds from a different molecule (%1) to "
            "a set of bonds from molecule %2.")
                .arg(bond.data().number())
                .arg(this->data().number()),
                    CODELOC);
    }

    auto atom0 = this->data().info().atomIdx(bond.ID().atom0());
    auto atom1 = this->data().info().atomIdx(bond.ID().atom1());

    if (atom0 > atom1)
        qSwap(atom0, atom1);

    if (atom0 == atom1)
        // cannot add bonds to the same atom
        return *this;

    SelectorBond ret(*this);

    ret.bnds.append(BondID(atom0, atom1));

    return ret;
}

MolViewPtr SelectorBond::operator[](int i) const
{
    return this->operator()(i);
}

MolViewPtr SelectorBond::operator[](const SireBase::Slice &slice) const
{
    return this->operator()(slice);
}

MolViewPtr SelectorBond::operator[](const QList<qint64> &idxs) const
{
    return this->operator()(idxs);
}

MolViewPtr SelectorBond::operator[](const BondID &bond) const
{
    return this->operator()(bond);
}

Bond SelectorBond::operator()(int i) const
{
    auto bnd = bnds.at(Index(i).map(bnds.count()));
    return Bond(this->data(), bnd);
}

SelectorBond SelectorBond::operator()(const SireBase::Slice &slice) const
{
    SelectorBond ret(*this);
    ret.bnds.clear();

    for (auto it = slice.begin(bnds.count()); not it.atEnd(); it.next())
    {
        ret.bnds.append(this->bnds.at(it.value()));
    }

    return ret;
}

SelectorBond SelectorBond::operator()(const QList<qint64> &idxs) const
{
    SelectorBond ret(*this);
    ret.bnds.clear();

    for (const auto &idx : idxs)
    {
        ret.bnds.append(this->bnds.at(Index(idx).map(this->bnds.count())));
    }

    return ret;
}

SelectorBond SelectorBond::operator()(int i, int j) const
{
    i = Index(i).map(bnds.count());
    j = Index(j).map(bnds.count());

    SelectorBond ret(*this);
    ret.bnds.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret.bnds.append(this->bnds.at(i));
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret.bnds.append(this->bnds.at(i));
        }
    }

    return ret;
}

SelectorBond SelectorBond::operator()(const BondID &bond) const
{
    auto atom0s = this->data().info().map(bond.atom0());
    auto atom1s = this->data().info().map(bond.atom1());

    SelectorBond ret(*this);
    ret.bnds.clear();

    for (const auto &atom0 : atom0s)
    {
        for (const auto &atom1 : atom1s)
        {
            auto a0 = atom0;
            auto a1 = atom1;

            if (a0 > a1)
                qSwap(a0, a1);

            BondID bond(a0, a1);

            for (const auto &b : bnds)
            {
                if (b == bond)
                    ret.bnds.append(b);
            }
        }
    }

    return ret;
}

MolViewPtr SelectorBond::toSelector() const
{
    return MolViewPtr(*this);
}

QList<MolViewPtr> SelectorBond::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(bnds.count());

    auto d = this->data();

    for (const auto &bnd : bnds)
    {
        l.append(MolViewPtr(new Bond(d, bnd)));
    }

    return l;
}

SelectorBond SelectorBond::add(const SelectorBond &other) const
{
    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;

    MoleculeView::assertSameMolecule(other);

    SelectorBond ret(*this);

    for (const auto &bond : other.bnds)
    {
        if (not this->bnds.contains(bond))
        {
            ret.bnds.append(bond);
        }
    }

    return ret;
}

SelectorBond SelectorBond::intersection(const SelectorBond &other) const
{
    if (this->isEmpty() or other.isEmpty())
    {
        return SelectorBond();
    }

    MoleculeView::assertSameMolecule(other);

    SelectorBond ret(*this);
    ret.bnds.clear();

    for (const auto &bond : this->bnds)
    {
        if (ret.bnds.contains(bond))
        {
            ret.bnds.append(bond);
        }
    }

    return ret;
}

SelectorBond SelectorBond::invert(const PropertyMap &map) const
{
    auto s = SelectorBond(this->molecule(), map);

    SelectorBond ret(*this);
    ret.bnds.clear();

    for (const auto &bond : s.bnds)
    {
        if (not this->bnds.contains(bond))
        {
            ret.bnds.append(bond);
        }
    }

    return ret;
}

SelectorBond SelectorBond::invert() const
{
    return this->invert(PropertyMap());
}

QList<BondID> SelectorBond::IDs() const
{
    return bnds;
}

bool SelectorBond::isEmpty() const
{
    return this->bnds.isEmpty();
}

bool SelectorBond::selectedAll() const
{
    return this->selection().selectedAll();
}

AtomSelection SelectorBond::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    for (const auto &bnd : bnds)
    {
        s = s.select(bnd.atom0());
        s = s.select(bnd.atom1());
    }

    return s;
}

bool SelectorBond::hasProperty(const PropertyName &key) const
{
    for (int i=0; i<this->count(); ++i)
    {
        if (this->operator()(i).hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorBond::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool SelectorBond::hasMetadata(const PropertyName &key,
                               const PropertyName &metakey) const
{
    return false;
}

QStringList SelectorBond::propertyKeys() const
{
    QSet<QString> keys;

    for (int i=0; i<this->count(); ++i)
    {
        for (const auto &k : this->operator()(i).propertyKeys())
        {
            keys.insert(k);
        }
    }

    //QStringList ret(keys.constBegin(), keys.constEnd());
    QStringList ret = keys.values();

    return ret;
}

QStringList SelectorBond::metadataKeys() const
{
    return QStringList();
}

QStringList SelectorBond::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

QList<Properties> SelectorBond::properties() const
{
    QList<Properties> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).properties());
    }

    return props;
}

Mover<SelectorBond> SelectorBond::move() const
{
    return Mover<SelectorBond>(*this);
}

Evaluator SelectorBond::evaluate() const
{
    return Evaluator(*this);
}

QList<PropertyPtr> SelectorBond::property(const PropertyName &key) const
{
    bool has_prop = false;

    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        try
        {
            props.append(this->operator()(i).property(key));
            has_prop = true;
        }
        catch(const SireError::exception&)
        {
            props.append(NullProperty());
        }
    }

    if (not has_prop)
        throw SireBase::missing_property(QObject::tr(
            "None of the bonds in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<PropertyPtr> SelectorBond::property(const PropertyName &key,
                                          const Property &default_value) const
{
    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).property(key, default_value));
    }

    return props;
}

QList<Length> SelectorBond::lengths(const PropertyMap &map) const
{
    QList<Length> l;

    for (int i=0; i<this->count(); ++i)
    {
        l.append(this->operator()(i).length(map));
    }

    return l;
}

QList<Length> SelectorBond::lengths() const
{
    return this->lengths(PropertyMap());
}

QList<Length> SelectorBond::measures() const
{
    return this->lengths();
}

QList<Length> SelectorBond::measures(const PropertyMap &map) const
{
    return this->lengths(map);
}

QList<Expression> SelectorBond::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<Expression> SelectorBond::potentials(const PropertyMap &map) const
{
    QList<Expression> p;

    for (int i=0; i<this->count(); ++i)
    {
        p.append(this->operator()(i).potential(map));
    }

    return p;
}

QList<MolarEnergy> SelectorBond::energies() const
{
    return this->energies(PropertyMap());
}

QList<MolarEnergy> SelectorBond::energies(const PropertyMap &map) const
{
    QList<MolarEnergy> nrgs;

    for (int i=0; i<this->count(); ++i)
    {
        nrgs.append(this->operator()(i).energy(map));
    }

    return nrgs;
}

MolarEnergy SelectorBond::energy() const
{
    return this->energy(PropertyMap());
}

MolarEnergy SelectorBond::energy(const PropertyMap &map) const
{
    MolarEnergy nrg(0);

    for (int i=0; i<this->count(); ++i)
    {
        nrg += this->operator()(i).energy(map);
    }

    return nrg;
}

namespace SireMol
{
    template class Mover<SelectorBond>;
}
