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

#include "selectorimproper.h"

#include "fouratomfunctions.h"

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

static const RegisterMetaType<SelectorImproper> r_simproper;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const SelectorImproper &simproper)
{
    writeHeader(ds, r_simproper, 1);
    SharedDataStream sds(ds);
    sds << simproper.imps << static_cast<const MoleculeView&>(simproper);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, SelectorImproper &simproper)
{
    auto v = readHeader(ds, r_simproper);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> simproper.imps >> static_cast<MoleculeView&>(simproper);
    }
    else
        throw SireStream::version_error(v, "1", r_simproper, CODELOC);

    return ds;
}

SelectorImproper::SelectorImproper() : ConcreteProperty<SelectorImproper, MoleculeView>()
{}

QList<ImproperID> _get_impropers(const MoleculeData &moldata,
                                 const PropertyMap &map)
{
    QList<ImproperID> impropers;

    if (moldata.hasProperty(map["improper"]))
    {
        auto funcs = moldata.property(map["improper"]).asA<FourAtomFunctions>();

        const auto &molinfo = moldata.info();

        for (const auto &func : funcs.potentials())
        {
            auto atomidx0 = molinfo.atomIdx(func.atom0());
            auto atomidx1 = molinfo.atomIdx(func.atom1());
            auto atomidx2 = molinfo.atomIdx(func.atom2());
            auto atomidx3 = molinfo.atomIdx(func.atom3());

            impropers.append(ImproperID(atomidx0, atomidx1,
                                        atomidx2, atomidx3));
        }
    }

    return impropers;
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto impropers = _get_impropers(mol.data(), map);

    if (mol.selectedAll())
    {
        imps = impropers;
    }
    else
    {
        const auto s = mol.selection();

        for (const auto &improper : impropers)
        {
            if (s.selected(improper.atom0()) and
                s.selected(improper.atom1()) and
                s.selected(improper.atom2()) and
                s.selected(improper.atom3()))
            {
                imps.append(improper);
            }
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const QList<ImproperID> &impropers)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    const auto s = mol.selection();

    for (const auto &improper : impropers)
    {
        ImproperID i(mol.data().info().atomIdx(improper.atom0()),
                     mol.data().info().atomIdx(improper.atom1()),
                     mol.data().info().atomIdx(improper.atom2()),
                     mol.data().info().atomIdx(improper.atom3()));

        if (s.selected(i.atom0()) and
            s.selected(i.atom1()) and
            s.selected(i.atom2()) and
            s.selected(i.atom3()))
        {
            imps.append(i);
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeData &moldata,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(moldata), map));
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const ImproperID &improper,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(improper.atom0());
    auto atoms1 = mol.data().info().map(improper.atom1());
    auto atoms2 = mol.data().info().map(improper.atom2());
    auto atoms3 = mol.data().info().map(improper.atom3());

    QList<ImproperID> all_impropers = _get_impropers(mol.data(), map);

    QSet<ImproperID> seen_impropers;

    QList<ImproperID> impropers;

    for (const auto &atom0 : atoms0)
    {
        for (const auto &atom1 : atoms1)
        {
            for (const auto &atom2 : atoms2)
            {
                for (const auto &atom3 : atoms3)
                {
                    auto atomidx0 = atom0;
                    auto atomidx1 = atom1;
                    auto atomidx2 = atom2;
                    auto atomidx3 = atom3;

                    ImproperID i(atomidx0, atomidx1, atomidx2, atomidx3);

                    if (atomidx0 != atomidx1 and
                        atomidx0 != atomidx2 and
                        atomidx0 != atomidx3 and
                        atomidx1 != atomidx2 and
                        atomidx1 != atomidx3 and
                        atomidx2 != atomidx3 and
                        not seen_impropers.contains(i))
                    {
                        // is this an improper?
                        for (const auto &imp : all_impropers)
                        {
                            if (imp == i)
                            {
                                seen_impropers.insert(i);
                                impropers.append(i);
                            }
                        }
                    }
                }
            }
        }

        if (mol.selectedAll())
        {
            imps = impropers;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &improper : impropers)
            {
                if (s.selected(improper.atom0()) and
                    s.selected(improper.atom1()) and
                    s.selected(improper.atom2()) and
                    s.selected(improper.atom3()))
                {
                    imps.append(improper);
                }
            }
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto all_impropers = _get_impropers(mol.data(), map);

    for (const auto &imp : all_impropers)
    {
        for (const auto &atomidx : mol.data().info().map(atom))
        {
            if (imp.atom0() == atomidx or
                imp.atom1() == atomidx or
                imp.atom2() == atomidx or
                imp.atom3() == atomidx)
            {
                imps.append(imp);
            }
        }
    }
}

bool _contains(const ImproperID &improper, const AtomIdx &atom)
{
    return atom == improper.atom0() or
           atom == improper.atom1() or
           atom == improper.atom2() or
           atom == improper.atom3();
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);

    auto impropers = _get_impropers(mol.data(), map);

    for (const auto &improper : impropers)
    {
        bool added = false;

        for (const auto &a0 : atoms0)
        {
            if (_contains(improper, a0))
            {
                for (const auto &a1 : atoms1)
                {
                    if (a0 != a1 and _contains(improper, a1))
                    {
                        imps.append(improper);
                        added = true;
                        break;
                    }
                }
            }

            if (added)
                break;
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);
    auto atoms2 = mol.data().info().map(atom2);

    auto impropers = _get_impropers(mol.data(), map);

    for (const auto &improper : impropers)
    {
        bool added = false;

        for (const auto &a0 : atoms0)
        {
            if (_contains(improper, a0))
            {
                for (const auto &a1 : atoms1)
                {
                    if (a0 != a1 and _contains(improper, a1))
                    {
                        for (const auto &a2 : atoms2)
                        {
                            if (a0 != a2 and a1 != a2 and _contains(improper, a2))
                            {
                                imps.append(improper);
                                added = true;
                                break;
                            }
                        }
                    }

                    if (added)
                        break;
                }
            }

            if (added)
                break;
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const AtomID &atom3,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(mol,
                                     ImproperID(atom0, atom1, atom2, atom3),
                                     map));
}

SelectorImproper::SelectorImproper(const MoleculeData &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(mol), atom, map));
}

SelectorImproper::SelectorImproper(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(mol), atom0, atom1, map));
}

SelectorImproper::SelectorImproper(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(mol), atom0, atom1, atom2, map));
}

SelectorImproper::SelectorImproper(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const AtomID &atom3,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(mol), atom0, atom1,
                                     atom2, atom3, map));
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(atoms)
{
    for (const auto &improper : _get_impropers(atoms.data(), map))
    {
        for (int i=0; i<atoms.count(); ++i)
        {
            if (_contains(improper, atoms(i).index()))
            {
                imps.append(improper);
                break;
            }
        }
    }
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a improper from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    for (const auto &improper : _get_impropers(atoms0.data(), map))
    {
        bool found = false;

        for (int i=0; i<atoms0.count(); ++i)
        {
            const auto atom0 = atoms0(i).index();

            if (_contains(improper, atom0))
            {
                for (int j=0; j<atoms1.count(); ++j)
                {
                    const auto atom1 = atoms1(j).index();

                    if (atom0 != atom1 and _contains(improper, atom1))
                    {
                        imps.append(improper);
                        found = true;
                        break;
                    }
                }
            }

            if (found)
                break;
        }
    }
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const PropertyMap &map)
              : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and atoms0.isSameMolecule(atoms2)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Improper from atoms in the same molecule. "
            "%1, %2 and %3 are from different molecules (%4, %5 and %6)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString()), CODELOC);

    for (const auto &improper : _get_impropers(atoms0.data(), map))
    {
        bool found = false;

        for (int i=0; i<atoms0.count(); ++i)
        {
            const auto atom0 = atoms0(i).index();

            if (_contains(improper, atom0))
            {
                for (int j=0; j<atoms1.count(); ++j)
                {
                    const auto atom1 = atoms1(j).index();

                    if (atom0 != atom1 and _contains(improper, atom1))
                    {
                        for (int k=0; k<atoms2.count(); ++k)
                        {
                            const auto atom2 = atoms2(k).index();

                            if (atom0 != atom2 and atom1 != atom2 and _contains(improper, atom2))
                            {
                                imps.append(improper);
                                found = true;
                                break;
                            }
                        }
                    }

                    if (found)
                        break;
                }
            }

            if (found)
                break;
        }
    }
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const Selector<Atom> &atoms3,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and
             atoms0.isSameMolecule(atoms2) and
             atoms0.isSameMolecule(atoms3)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an Improper from atoms in the same molecule. "
            "%1, %2, %3 and %4 are from different molecules (%5, %6, %7 and %8)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString()).arg(atoms3.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString())
                .arg(atoms3.molecule().toString()), CODELOC);

    for (const auto &improper : _get_impropers(atoms0.data(), map))
    {
        bool found = false;

        for (int i=0; i<atoms0.count(); ++i)
        {
            const auto atom0 = atoms0(i).index();

            if (_contains(improper, atom0))
            {
                for (int j=0; j<atoms1.count(); ++j)
                {
                    const auto atom1 = atoms1(j).index();

                    if (atom0 != atom1 and _contains(improper, atom1))
                    {
                        for (int k=0; k<atoms2.count(); ++k)
                        {
                            const auto atom2 = atoms2(k).index();

                            if (atom0 != atom2 and atom1 != atom2 and _contains(improper, atom2))
                            {
                                for (int l=0; l<atoms3.count(); ++l)
                                {
                                    const auto atom3 = atoms3(l).index();

                                    if (atom0 != atom3 and atom1 != atom3 and atom2 != atom3 and _contains(improper, atom3))
                                    {
                                        imps.append(improper);
                                        found = true;
                                        break;
                                    }
                                }
                            }

                            if (found)
                                break;
                        }
                    }

                    if (found)
                        break;
                }
            }

            if (found)
                break;
        }
    }
}

SelectorImproper::SelectorImproper(const SelectorImproper &other)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(other),
                   imps(other.imps)
{}

SelectorImproper::~SelectorImproper()
{}

const char* SelectorImproper::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorImproper>());
}

SelectorImproper& SelectorImproper::operator=(const SelectorImproper &other)
{
    if (this != &other)
    {
        imps = other.imps;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorImproper::operator==(const SelectorImproper &other) const
{
    return imps == other.imps and MoleculeView::operator==(other);
}

bool SelectorImproper::operator!=(const SelectorImproper &other) const
{
    return not operator==(other);
}

int SelectorImproper::count() const
{
    return this->imps.count();
}

int SelectorImproper::size() const
{
    return this->count();
}

int SelectorImproper::nViews() const
{
    return this->count();
}

QString SelectorImproper::toString() const
{
    if (this->isNull())
        return QObject::tr("SelectorImproper::null");
    else if (this->isEmpty())
        return QObject::tr("SelectorImproper::empty");

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

    return QObject::tr("SelectorImproper( size=%1\n%2\n)")
                .arg(this->count()).arg(parts.join("\n"));
}

SelectorImproper SelectorImproper::add(const Improper &improper) const
{
    if (improper.isNull())
        return *this;

    if (this->isEmpty())
    {
        return SelectorImproper(improper);
    }

    if (improper.data().number() != this->data().number())
    {
        throw SireError::incompatible_error(QObject::tr(
            "You cannot add Impropers from a different molecule (%1) to "
            "a set of Impropers from molecule %2.")
                .arg(improper.data().number())
                .arg(this->data().number()),
                    CODELOC);
    }

    auto atom0 = this->data().info().atomIdx(improper.ID().atom0());
    auto atom1 = this->data().info().atomIdx(improper.ID().atom1());
    auto atom2 = this->data().info().atomIdx(improper.ID().atom2());
    auto atom3 = this->data().info().atomIdx(improper.ID().atom3());


    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
        atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        // cannot add Impropers to the same atom
        return *this;

    SelectorImproper ret(*this);

    ret.imps.append(ImproperID(atom0, atom1, atom2, atom3));

    return ret;
}

MolViewPtr SelectorImproper::operator[](int i) const
{
    return this->operator()(i);
}

MolViewPtr SelectorImproper::operator[](const SireBase::Slice &slice) const
{
    return this->operator()(slice);
}

MolViewPtr SelectorImproper::operator[](const QList<qint64> &idxs) const
{
    return this->operator()(idxs);
}

MolViewPtr SelectorImproper::operator[](const ImproperID &Improper) const
{
    return this->operator()(Improper);
}

SireMM::Improper SelectorImproper::operator()(int i) const
{
    auto imp = imps.at(Index(i).map(imps.count()));
    return SireMM::Improper(this->data(), imp);
}

SelectorImproper SelectorImproper::operator()(const SireBase::Slice &slice) const
{
    SelectorImproper ret(*this);
    ret.imps.clear();

    for (auto it = slice.begin(imps.count()); not it.atEnd(); it.next())
    {
        ret.imps.append(this->imps.at(it.value()));
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(const QList<qint64> &idxs) const
{
    SelectorImproper ret(*this);
    ret.imps.clear();

    for (const auto &idx : idxs)
    {
        ret.imps.append(this->imps.at(Index(idx).map(this->imps.count())));
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(int i, int j) const
{
    i = Index(i).map(imps.count());
    j = Index(j).map(imps.count());

    SelectorImproper ret(*this);
    ret.imps.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret.imps.append(this->imps.at(i));
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret.imps.append(this->imps.at(i));
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(const ImproperID &improper) const
{
    auto atom0s = this->data().info().map(improper.atom0());
    auto atom1s = this->data().info().map(improper.atom1());
    auto atom2s = this->data().info().map(improper.atom2());
    auto atom3s = this->data().info().map(improper.atom3());

    SelectorImproper ret(*this);
    ret.imps.clear();

    for (const auto &imp : imps)
    {
        for (const auto &atom0 : atom0s)
        {
            for (const auto &atom1 : atom1s)
            {
                for (const auto &atom2 : atom2s)
                {
                    for (const auto &atom3 : atom3s)
                    {
                        auto a0 = atom0;
                        auto a1 = atom1;
                        auto a2 = atom2;
                        auto a3 = atom3;

                        ImproperID improper(a0, a1, a2, a3);

                        if (imp == improper)
                            ret.imps.append(imp);
                    }
                }
            }
        }
    }

    return ret;
}

MolViewPtr SelectorImproper::toSelector() const
{
    return MolViewPtr(*this);
}

QList<MolViewPtr> SelectorImproper::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(imps.count());

    auto d = this->data();

    for (const auto &imp : imps)
    {
        l.append(MolViewPtr(new Improper(d, imp)));
    }

    return l;
}

SelectorImproper SelectorImproper::add(const SelectorImproper &other) const
{
    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;

    MoleculeView::assertSameMolecule(other);

    SelectorImproper ret(*this);

    for (const auto &improper : other.imps)
    {
        if (not this->imps.contains(improper))
        {
            ret.imps.append(improper);
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::intersection(const SelectorImproper &other) const
{
    if (this->isEmpty() or other.isEmpty())
    {
        return SelectorImproper();
    }

    MoleculeView::assertSameMolecule(other);

    SelectorImproper ret(*this);
    ret.imps.clear();

    for (const auto &improper : this->imps)
    {
        if (ret.imps.contains(improper))
        {
            ret.imps.append(improper);
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::invert(const PropertyMap &map) const
{
    auto s = SelectorImproper(this->molecule(), map);

    SelectorImproper ret(*this);
    ret.imps.clear();

    for (const auto &improper : s.imps)
    {
        if (not this->imps.contains(improper))
        {
            ret.imps.append(improper);
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::invert() const
{
    return this->invert(PropertyMap());
}

QList<ImproperID> SelectorImproper::IDs() const
{
    return imps;
}

bool SelectorImproper::isEmpty() const
{
    return this->imps.isEmpty();
}

bool SelectorImproper::selectedAll() const
{
    return this->selection().selectedAll();
}

AtomSelection SelectorImproper::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    for (const auto &imp : imps)
    {
        s = s.select(imp.atom0());
        s = s.select(imp.atom1());
        s = s.select(imp.atom2());
        s = s.select(imp.atom3());
    }

    return s;
}

bool SelectorImproper::hasProperty(const PropertyName &key) const
{
    for (int i=0; i<this->count(); ++i)
    {
        if (this->operator()(i).hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorImproper::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool SelectorImproper::hasMetadata(const PropertyName &key,
                                   const PropertyName &metakey) const
{
    return false;
}

QStringList SelectorImproper::propertyKeys() const
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

QStringList SelectorImproper::metadataKeys() const
{
    return QStringList();
}

QStringList SelectorImproper::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

QList<Properties> SelectorImproper::properties() const
{
    QList<Properties> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).properties());
    }

    return props;
}

Mover<SelectorImproper> SelectorImproper::move() const
{
    return Mover<SelectorImproper>(*this);
}

Evaluator SelectorImproper::evaluate() const
{
    return Evaluator(*this);
}

QList<PropertyPtr> SelectorImproper::property(const PropertyName &key) const
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
        catch(SireError::exception&)
        {
            props.append(NullProperty());
        }
    }

    if (not has_prop)
        throw SireBase::missing_property(QObject::tr(
            "None of the impropers in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<PropertyPtr> SelectorImproper::property(const PropertyName &key,
                                              const Property &default_value) const
{
    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).property(key, default_value));
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorImproper::sizes(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> s;

    for (int i=0; i<this->count(); ++i)
    {
        s.append(this->operator()(i).size(map));
    }

    return s;
}

QList<SireUnits::Dimension::Angle> SelectorImproper::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<Expression> SelectorImproper::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<Expression> SelectorImproper::potentials(const PropertyMap &map) const
{
    QList<Expression> p;

    for (int i=0; i<this->count(); ++i)
    {
        p.append(this->operator()(i).potential(map));
    }

    return p;
}

QList<SireUnits::Dimension::MolarEnergy> SelectorImproper::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::MolarEnergy> SelectorImproper::energies(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::MolarEnergy> nrgs;

    for (int i=0; i<this->count(); ++i)
    {
        nrgs.append(this->operator()(i).energy(map));
    }

    return nrgs;
}

SireUnits::Dimension::MolarEnergy SelectorImproper::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy SelectorImproper::energy(const PropertyMap &map) const
{
    SireUnits::Dimension::MolarEnergy nrg(0);

    for (int i=0; i<this->count(); ++i)
    {
        nrg += this->operator()(i).energy(map);
    }

    return nrg;
}

namespace SireMol
{
    template class Mover<SelectorImproper>;
}
