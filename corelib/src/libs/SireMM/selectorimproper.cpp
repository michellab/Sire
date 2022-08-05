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
    QList<ImproperID> imps;

    if (mol.data().hasProperty(map["improper"]))
    {
        auto funcs = mol.data().property(map["improper"]).asA<FourAtomFunctions>();

        const auto &molinfo = mol.data().info();

        QList<ImproperID> impropers;

        for (const auto func : funcs.potentials())
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
                        not seen_Impropers.contains(i))
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
                if (s.selected(Improper.atom0()) and
                    s.selected(Improper.atom1()) and
                    s.selected(Improper.atom2()) and
                    s.selected(Improper.atom3()))
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
        for (const auto &atomidx : mol.data().info().atomIdxs(atom))
        {
            if (imp.atom0() == atomidx or
                imp.atom1() == atomidx or
                imp.atom2() == atomidx or
                imp.atom3() == atomidx)
            {
                ims.append(imp);
            }
        }
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
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
    if (atoms.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<ImproperID> seen_dihs;

        QList<ImproperID> Impropers;

        for (int i=0; i<atoms.count(); ++i)
        {
            for (const auto &d : c.getImpropers(atoms(i).index()))
            {
                auto atomidx0 = atoms.data().info().atomIdx(d.atom0());
                auto atomidx1 = atoms.data().info().atomIdx(d.atom1());
                auto atomidx2 = atoms.data().info().atomIdx(d.atom2());
                auto atomidx3 = atoms.data().info().atomIdx(d.atom3());

                if (atomidx0 > atomidx3)
                {
                    qSwap(atomidx0, atomidx3);
                    qSwap(atomidx1, atomidx2);
                }

                ImproperID dih(atomidx0, atomidx1, atomidx2, atomidx3);

                if (atomidx0 != atomidx1 and
                    atomidx0 != atomidx2 and
                    atomidx0 != atomidx3 and
                    atomidx1 != atomidx2 and
                    atomidx1 != atomidx3 and
                    atomidx2 != atomidx3 and
                    not seen_dihs.contains(dih))
                {
                    seen_dihs.insert(dih);
                    Impropers.append(dih);
                }
            }
        }

        dihs = Impropers;
    }
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const PropertyMap &map)
              : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Improper from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<ImproperID> seen_dihs;

        QList<ImproperID> Impropers;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                auto atomidx0 = atoms0(i).index();
                auto atomidx1 = atoms1(j).index();

                auto dihs = c.getImpropers(atomidx0, atomidx1);

                for (const auto &dih : dihs)
                {
                    if (not seen_dihs.contains(dih))
                    {
                        seen_dihs.insert(dih);
                        Impropers.append(dih);
                    }
                }
            }
        }

        dihs = Impropers;
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

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<ImproperID> seen_dihs;

        QList<ImproperID> Impropers;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                for (int k=0; k<atoms2.count(); ++k)
                {
                    auto atomidx0 = atoms0(i).index();
                    auto atomidx1 = atoms1(j).index();
                    auto atomidx2 = atoms2(k).index();

                    auto dihs = c.getImpropers(atomidx0, atomidx1, atomidx2);

                    for (const auto &dih : dihs)
                    {
                        if (not seen_dihs.contains(dih))
                        {
                            seen_dihs.insert(dih);
                            Impropers.append(dih);
                        }
                    }
                }
            }
        }

        dihs = Impropers;
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

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<ImproperID> seen_dihs;

        QList<ImproperID> Impropers;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                for (int k=0; k<atoms2.count(); ++k)
                {
                    for (int l=0; l<atoms3.count(); ++l)
                    {
                        auto atomidx0 = atoms0(i).index();
                        auto atomidx1 = atoms1(j).index();
                        auto atomidx2 = atoms2(k).index();
                        auto atomidx3 = atoms3(l).index();

                        if (c.areImpropered(atomidx0, atomidx1, atomidx2, atomidx3))
                        {
                            if (atomidx0 > atomidx3)
                            {
                                qSwap(atomidx0, atomidx3);
                                qSwap(atomidx1, atomidx2);
                            }

                            ImproperID d(atomidx0, atomidx1, atomidx2, atomidx3);

                            if (not seen_dihs.contains(d))
                            {
                                seen_dihs.insert(d);
                                Impropers.append(d);
                            }
                        }
                    }
                }
            }
        }

        dihs = Impropers;
    }
}

SelectorImproper::SelectorImproper(const SelectorImproper &other)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(other),
                   dihs(other.dihs)
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
        dihs = other.dihs;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorImproper::operator==(const SelectorImproper &other) const
{
    return dihs == other.dihs and MoleculeView::operator==(other);
}

bool SelectorImproper::operator!=(const SelectorImproper &other) const
{
    return not operator==(other);
}

int SelectorImproper::count() const
{
    return this->dihs.count();
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

SelectorImproper SelectorImproper::add(const Improper &Improper) const
{
    if (Improper.isNull())
        return *this;

    if (this->isEmpty())
    {
        return SelectorImproper(Improper);
    }

    if (Improper.data().number() != this->data().number())
    {
        throw SireError::incompatible_error(QObject::tr(
            "You cannot add Impropers from a different molecule (%1) to "
            "a set of Impropers from molecule %2.")
                .arg(Improper.data().number())
                .arg(this->data().number()),
                    CODELOC);
    }

    auto atom0 = this->data().info().atomIdx(Improper.ID().atom0());
    auto atom1 = this->data().info().atomIdx(Improper.ID().atom1());
    auto atom2 = this->data().info().atomIdx(Improper.ID().atom2());
    auto atom3 = this->data().info().atomIdx(Improper.ID().atom3());

    if (atom0 > atom3)
    {
        qSwap(atom0, atom3);
        qSwap(atom1, atom2);
    }

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
        atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        // cannot add Impropers to the same atom
        return *this;

    SelectorImproper ret(*this);

    ret.dihs.append(ImproperID(atom0, atom1, atom2, atom3));

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
    auto dih = dihs.at(Index(i).map(dihs.count()));
    return SireMM::Improper(this->data(), dih);
}

SelectorImproper SelectorImproper::operator()(const SireBase::Slice &slice) const
{
    SelectorImproper ret(*this);
    ret.dihs.clear();

    for (auto it = slice.begin(dihs.count()); not it.atEnd(); it.next())
    {
        ret.dihs.append(this->dihs.at(it.value()));
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(const QList<qint64> &idxs) const
{
    SelectorImproper ret(*this);
    ret.dihs.clear();

    for (const auto &idx : idxs)
    {
        ret.dihs.append(this->dihs.at(Index(idx).map(this->dihs.count())));
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(int i, int j) const
{
    i = Index(i).map(dihs.count());
    j = Index(j).map(dihs.count());

    SelectorImproper ret(*this);
    ret.dihs.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret.dihs.append(this->dihs.at(i));
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret.dihs.append(this->dihs.at(i));
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::operator()(const ImproperID &Improper) const
{
    auto atom0s = this->data().info().map(Improper.atom0());
    auto atom1s = this->data().info().map(Improper.atom1());
    auto atom2s = this->data().info().map(Improper.atom2());
    auto atom3s = this->data().info().map(Improper.atom3());

    SelectorImproper ret(*this);
    ret.dihs.clear();

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

                    if (a0 > a3)
                    {
                        qSwap(a0, a3);
                        qSwap(a1, a2);
                    }

                    ImproperID Improper(a0, a1, a2, a3);

                    for (const auto &d : dihs)
                    {
                        if (d == Improper)
                            ret.dihs.append(d);
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
    l.reserve(dihs.count());

    auto d = this->data();

    for (const auto &dih : dihs)
    {
        l.append(MolViewPtr(new Improper(d, dih)));
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

    for (const auto &Improper : other.dihs)
    {
        if (not this->dihs.contains(Improper))
        {
            ret.dihs.append(Improper);
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
    ret.dihs.clear();

    for (const auto &Improper : this->dihs)
    {
        if (ret.dihs.contains(Improper))
        {
            ret.dihs.append(Improper);
        }
    }

    return ret;
}

SelectorImproper SelectorImproper::invert(const PropertyMap &map) const
{
    auto s = SelectorImproper(this->molecule(), map);

    SelectorImproper ret(*this);
    ret.dihs.clear();

    for (const auto &Improper : s.dihs)
    {
        if (not this->dihs.contains(Improper))
        {
            ret.dihs.append(Improper);
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
    return dihs;
}

bool SelectorImproper::isEmpty() const
{
    return this->dihs.isEmpty();
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

    for (const auto &dih : dihs)
    {
        s = s.select(dih.atom0());
        s = s.select(dih.atom1());
        s = s.select(dih.atom2());
        s = s.select(dih.atom3());
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
            "None of the Impropers in this container have a property called %1.")
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
