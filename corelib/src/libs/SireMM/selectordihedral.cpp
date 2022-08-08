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

#include "selectordihedral.h"

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

static const RegisterMetaType<SelectorDihedral> r_sdihedral;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const SelectorDihedral &sdihedral)
{
    writeHeader(ds, r_sdihedral, 1);
    SharedDataStream sds(ds);
    sds << sdihedral.dihs << static_cast<const MoleculeView&>(sdihedral);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, SelectorDihedral &sdihedral)
{
    auto v = readHeader(ds, r_sdihedral);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> sdihedral.dihs >> static_cast<MoleculeView&>(sdihedral);
    }
    else
        throw SireStream::version_error(v, "1", r_sdihedral, CODELOC);

    return ds;
}

SelectorDihedral::SelectorDihedral() : ConcreteProperty<SelectorDihedral, MoleculeView>()
{}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto dihedrals = c.getDihedrals();

        if (mol.selectedAll())
        {
            dihs = dihedrals;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &dihedral : dihedrals)
            {
                if (s.selected(dihedral.atom0()) and
                    s.selected(dihedral.atom1()) and
                    s.selected(dihedral.atom2()) and
                    s.selected(dihedral.atom3()))
                {
                    dihs.append(dihedral);
                }
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const QList<DihedralID> &dihedrals)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    const auto s = mol.selection();

    for (const auto &dihedral : dihedrals)
    {
        DihedralID d(mol.data().info().atomIdx(dihedral.atom0()),
                     mol.data().info().atomIdx(dihedral.atom1()),
                     mol.data().info().atomIdx(dihedral.atom2()),
                     mol.data().info().atomIdx(dihedral.atom3()));

        if (s.selected(d.atom0()) and
            s.selected(d.atom1()) and
            s.selected(d.atom2()) and
            s.selected(d.atom3()))
        {
            dihs.append(d);
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeData &moldata,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(moldata), map));
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const DihedralID &dihedral,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(dihedral.atom0());
    auto atoms1 = mol.data().info().map(dihedral.atom1());
    auto atoms2 = mol.data().info().map(dihedral.atom2());
    auto atoms3 = mol.data().info().map(dihedral.atom3());

    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<DihedralID> seen_dihedrals;

        QList<DihedralID> dihedrals;

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

                        if (atomidx0 > atomidx3)
                        {
                            qSwap(atomidx0, atomidx3);
                            qSwap(atomidx1, atomidx2);
                        }

                        DihedralID d(atomidx0, atomidx1, atomidx2, atomidx3);

                        if (atomidx0 != atomidx1 and
                            atomidx1 != atomidx2 and
                            atomidx2 != atomidx3 and
                            c.areConnected(atomidx0, atomidx1) and
                            c.areConnected(atomidx1, atomidx2) and
                            c.areConnected(atomidx2, atomidx3) and
                            not seen_dihedrals.contains(d))
                        {
                            seen_dihedrals.insert(d);
                            dihedrals.append(d);
                        }
                    }
                }
            }
        }

        if (mol.selectedAll())
        {
            dihs = dihedrals;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &dihedral : dihedrals)
            {
                if (s.selected(dihedral.atom0()) and
                    s.selected(dihedral.atom1()) and
                    s.selected(dihedral.atom2()) and
                    s.selected(dihedral.atom3()))
                {
                    dihs.append(dihedral);
                }
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto dihedrals = c.getDihedrals(atom);

        if (mol.selectedAll())
        {
            dihs = dihedrals;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &dihedral : dihedrals)
            {
                if (s.selected(dihedral.atom0()) and
                    s.selected(dihedral.atom1()) and
                    s.selected(dihedral.atom2()) and
                    s.selected(dihedral.atom3()))
                {
                    dihs.append(dihedral);
                }
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto dihedrals = c.getDihedrals(atom0, atom1);

        if (mol.selectedAll())
        {
            dihs = dihedrals;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &dihedral : dihedrals)
            {
                if (s.selected(dihedral.atom0()) and
                    s.selected(dihedral.atom1()) and
                    s.selected(dihedral.atom2()) and
                    s.selected(dihedral.atom3()))
                {
                    dihs.append(dihedral);
                }
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto dihedrals = c.getDihedrals(atom0, atom1, atom2);

        if (mol.selectedAll())
        {
            dihs = dihedrals;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &dihedral : dihedrals)
            {
                if (s.selected(dihedral.atom0()) and
                    s.selected(dihedral.atom1()) and
                    s.selected(dihedral.atom2()) and
                    s.selected(dihedral.atom3()))
                {
                    dihs.append(dihedral);
                }
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const AtomID &atom3,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(mol,
                                     DihedralID(atom0, atom1, atom2, atom3),
                                     map));
}

SelectorDihedral::SelectorDihedral(const MoleculeData &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(mol), atom, map));
}

SelectorDihedral::SelectorDihedral(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(mol), atom0, atom1, map));
}

SelectorDihedral::SelectorDihedral(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(mol), atom0, atom1, atom2, map));
}

SelectorDihedral::SelectorDihedral(const MoleculeData &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const AtomID &atom3,
                                   const PropertyMap &map)
             : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(mol), atom0, atom1,
                                     atom2, atom3, map));
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms)
{
    if (atoms.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<DihedralID> seen_dihs;

        QList<DihedralID> dihedrals;

        for (int i=0; i<atoms.count(); ++i)
        {
            for (const auto &d : c.getDihedrals(atoms(i).index()))
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

                DihedralID dih(atomidx0, atomidx1, atomidx2, atomidx3);

                if (atomidx0 != atomidx1 and
                    atomidx1 != atomidx2 and
                    atomidx2 != atomidx3 and
                    not seen_dihs.contains(dih))
                {
                    seen_dihs.insert(dih);
                    dihedrals.append(dih);
                }
            }
        }

        dihs = dihedrals;
    }
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const PropertyMap &map)
              : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Dihedral from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<DihedralID> seen_dihs;

        QList<DihedralID> dihedrals;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                auto atomidx0 = atoms0(i).index();
                auto atomidx1 = atoms1(j).index();

                auto dihs = c.getDihedrals(atomidx0, atomidx1);

                for (const auto &dih : dihs)
                {
                    if (not seen_dihs.contains(dih))
                    {
                        seen_dihs.insert(dih);
                        dihedrals.append(dih);
                    }
                }
            }
        }

        dihs = dihedrals;
    }
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const PropertyMap &map)
              : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and atoms0.isSameMolecule(atoms2)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Dihedral from atoms in the same molecule. "
            "%1, %2 and %3 are from different molecules (%4, %5 and %6)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString()), CODELOC);

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<DihedralID> seen_dihs;

        QList<DihedralID> dihedrals;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                for (int k=0; k<atoms2.count(); ++k)
                {
                    auto atomidx0 = atoms0(i).index();
                    auto atomidx1 = atoms1(j).index();
                    auto atomidx2 = atoms2(k).index();

                    auto dihs = c.getDihedrals(atomidx0, atomidx1, atomidx2);

                    for (const auto &dih : dihs)
                    {
                        if (not seen_dihs.contains(dih))
                        {
                            seen_dihs.insert(dih);
                            dihedrals.append(dih);
                        }
                    }
                }
            }
        }

        dihs = dihedrals;
    }
}

bool _contains(const DihedralID &dih, const AtomIdx &atom)
{
    return dih.atom0() == atom or
           dih.atom1() == atom or
           dih.atom2() == atom or
           dih.atom3() == atom;
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const Selector<Atom> &atoms3,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and
             atoms0.isSameMolecule(atoms2) and
             atoms0.isSameMolecule(atoms3)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an Dihedral from atoms in the same molecule. "
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

        for (const auto &dihedral : c.getDihedrals())
        {
            bool found = false;

            for (int i=0; i<atoms0.count(); ++i)
            {
                auto atomidx0 = atoms0(i).index();

                if (not _contains(dihedral, atomidx0))
                    break;

                for (int j=0; j<atoms1.count(); ++j)
                {
                    auto atomidx1 = atoms1(j).index();

                    if (atomidx1 == atomidx0 or not _contains(dihedral, atomidx1))
                        break;

                    for (int k=0; k<atoms2.count(); ++k)
                    {
                        auto atomidx2 = atoms2(k).index();

                        if (atomidx2 == atomidx1 or not _contains(dihedral, atomidx2))
                            break;


                        for (int l=0; l<atoms3.count(); ++l)
                        {
                            auto atomidx3 = atoms3(l).index();

                            if (atomidx3 == atomidx2 or not _contains(dihedral, atomidx3))
                                break;

                            if (atomidx0 > atomidx3)
                            {
                                qSwap(atomidx0, atomidx3);
                                qSwap(atomidx1, atomidx2);
                            }

                            dihs.append(DihedralID(atomidx0, atomidx1, atomidx2, atomidx3));
                            found = true;
                            break;
                        }

                        if (found)
                            break;
                    }

                    if (found)
                        break;
                }

                if (found)
                    break;
            }
        }
    }
}

SelectorDihedral::SelectorDihedral(const SelectorDihedral &other)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(other),
                   dihs(other.dihs)
{}

SelectorDihedral::~SelectorDihedral()
{}

const char* SelectorDihedral::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorDihedral>());
}

SelectorDihedral& SelectorDihedral::operator=(const SelectorDihedral &other)
{
    if (this != &other)
    {
        dihs = other.dihs;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorDihedral::operator==(const SelectorDihedral &other) const
{
    return dihs == other.dihs and MoleculeView::operator==(other);
}

bool SelectorDihedral::operator!=(const SelectorDihedral &other) const
{
    return not operator==(other);
}

int SelectorDihedral::count() const
{
    return this->dihs.count();
}

int SelectorDihedral::size() const
{
    return this->count();
}

int SelectorDihedral::nViews() const
{
    return this->count();
}

QString SelectorDihedral::toString() const
{
    if (this->isNull())
        return QObject::tr("SelectorDihedral::null");
    else if (this->isEmpty())
        return QObject::tr("SelectorDihedral::empty");

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

    return QObject::tr("SelectorDihedral( size=%1\n%2\n)")
                .arg(this->count()).arg(parts.join("\n"));
}

SelectorDihedral SelectorDihedral::add(const Dihedral &dihedral) const
{
    if (dihedral.isNull())
        return *this;

    if (this->isEmpty())
    {
        return SelectorDihedral(dihedral);
    }

    if (dihedral.data().number() != this->data().number())
    {
        throw SireError::incompatible_error(QObject::tr(
            "You cannot add Dihedrals from a different molecule (%1) to "
            "a set of Dihedrals from molecule %2.")
                .arg(dihedral.data().number())
                .arg(this->data().number()),
                    CODELOC);
    }

    auto atom0 = this->data().info().atomIdx(dihedral.ID().atom0());
    auto atom1 = this->data().info().atomIdx(dihedral.ID().atom1());
    auto atom2 = this->data().info().atomIdx(dihedral.ID().atom2());
    auto atom3 = this->data().info().atomIdx(dihedral.ID().atom3());

    if (atom0 > atom3)
    {
        qSwap(atom0, atom3);
        qSwap(atom1, atom2);
    }

    if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
        atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        // cannot add Dihedrals to the same atom
        return *this;

    SelectorDihedral ret(*this);

    ret.dihs.append(DihedralID(atom0, atom1, atom2, atom3));

    return ret;
}

MolViewPtr SelectorDihedral::operator[](int i) const
{
    return this->operator()(i);
}

MolViewPtr SelectorDihedral::operator[](const SireBase::Slice &slice) const
{
    return this->operator()(slice);
}

MolViewPtr SelectorDihedral::operator[](const QList<qint64> &idxs) const
{
    return this->operator()(idxs);
}

MolViewPtr SelectorDihedral::operator[](const DihedralID &Dihedral) const
{
    return this->operator()(Dihedral);
}

SireMM::Dihedral SelectorDihedral::operator()(int i) const
{
    auto dih = dihs.at(Index(i).map(dihs.count()));
    return SireMM::Dihedral(this->data(), dih);
}

SelectorDihedral SelectorDihedral::operator()(const SireBase::Slice &slice) const
{
    SelectorDihedral ret(*this);
    ret.dihs.clear();

    for (auto it = slice.begin(dihs.count()); not it.atEnd(); it.next())
    {
        ret.dihs.append(this->dihs.at(it.value()));
    }

    return ret;
}

SelectorDihedral SelectorDihedral::operator()(const QList<qint64> &idxs) const
{
    SelectorDihedral ret(*this);
    ret.dihs.clear();

    for (const auto &idx : idxs)
    {
        ret.dihs.append(this->dihs.at(Index(idx).map(this->dihs.count())));
    }

    return ret;
}

SelectorDihedral SelectorDihedral::operator()(int i, int j) const
{
    i = Index(i).map(dihs.count());
    j = Index(j).map(dihs.count());

    SelectorDihedral ret(*this);
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

SelectorDihedral SelectorDihedral::operator()(const DihedralID &dihedral) const
{
    auto atom0s = this->data().info().map(dihedral.atom0());
    auto atom1s = this->data().info().map(dihedral.atom1());
    auto atom2s = this->data().info().map(dihedral.atom2());
    auto atom3s = this->data().info().map(dihedral.atom3());

    SelectorDihedral ret(*this);
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

                    DihedralID dihedral(a0, a1, a2, a3);

                    for (const auto &d : dihs)
                    {
                        if (d == dihedral)
                            ret.dihs.append(d);
                    }
                }
            }
        }
    }

    return ret;
}

MolViewPtr SelectorDihedral::toSelector() const
{
    return MolViewPtr(*this);
}

QList<MolViewPtr> SelectorDihedral::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(dihs.count());

    auto d = this->data();

    for (const auto &dih : dihs)
    {
        l.append(MolViewPtr(new Dihedral(d, dih)));
    }

    return l;
}

SelectorDihedral SelectorDihedral::add(const SelectorDihedral &other) const
{
    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;

    MoleculeView::assertSameMolecule(other);

    SelectorDihedral ret(*this);

    for (const auto &dihedral : other.dihs)
    {
        if (not this->dihs.contains(dihedral))
        {
            ret.dihs.append(dihedral);
        }
    }

    return ret;
}

SelectorDihedral SelectorDihedral::intersection(const SelectorDihedral &other) const
{
    if (this->isEmpty() or other.isEmpty())
    {
        return SelectorDihedral();
    }

    MoleculeView::assertSameMolecule(other);

    SelectorDihedral ret(*this);
    ret.dihs.clear();

    for (const auto &dihedral : this->dihs)
    {
        if (ret.dihs.contains(dihedral))
        {
            ret.dihs.append(dihedral);
        }
    }

    return ret;
}

SelectorDihedral SelectorDihedral::invert(const PropertyMap &map) const
{
    auto s = SelectorDihedral(this->molecule(), map);

    SelectorDihedral ret(*this);
    ret.dihs.clear();

    for (const auto &dihedral : s.dihs)
    {
        if (not this->dihs.contains(dihedral))
        {
            ret.dihs.append(dihedral);
        }
    }

    return ret;
}

SelectorDihedral SelectorDihedral::invert() const
{
    return this->invert(PropertyMap());
}

QList<DihedralID> SelectorDihedral::IDs() const
{
    return dihs;
}

bool SelectorDihedral::isEmpty() const
{
    return this->dihs.isEmpty();
}

bool SelectorDihedral::selectedAll() const
{
    return this->selection().selectedAll();
}

AtomSelection SelectorDihedral::selection() const
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

bool SelectorDihedral::hasProperty(const PropertyName &key) const
{
    for (int i=0; i<this->count(); ++i)
    {
        if (this->operator()(i).hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorDihedral::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool SelectorDihedral::hasMetadata(const PropertyName &key,
                                   const PropertyName &metakey) const
{
    return false;
}

QStringList SelectorDihedral::propertyKeys() const
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

QStringList SelectorDihedral::metadataKeys() const
{
    return QStringList();
}

QStringList SelectorDihedral::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

QList<Properties> SelectorDihedral::properties() const
{
    QList<Properties> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).properties());
    }

    return props;
}

Mover<SelectorDihedral> SelectorDihedral::move() const
{
    return Mover<SelectorDihedral>(*this);
}

Evaluator SelectorDihedral::evaluate() const
{
    return Evaluator(*this);
}

QList<PropertyPtr> SelectorDihedral::property(const PropertyName &key) const
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
            "None of the dihedrals in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<PropertyPtr> SelectorDihedral::property(const PropertyName &key,
                                              const Property &default_value) const
{
    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).property(key, default_value));
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorDihedral::sizes(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> s;

    for (int i=0; i<this->count(); ++i)
    {
        s.append(this->operator()(i).size(map));
    }

    return s;
}

QList<SireUnits::Dimension::Angle> SelectorDihedral::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<SireUnits::Dimension::Angle> SelectorDihedral::measures(const PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireUnits::Dimension::Angle> SelectorDihedral::measures() const
{
    return this->sizes();
}

QList<Expression> SelectorDihedral::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<Expression> SelectorDihedral::potentials(const PropertyMap &map) const
{
    QList<Expression> p;

    for (int i=0; i<this->count(); ++i)
    {
        p.append(this->operator()(i).potential(map));
    }

    return p;
}

QList<SireUnits::Dimension::MolarEnergy> SelectorDihedral::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::MolarEnergy> SelectorDihedral::energies(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::MolarEnergy> nrgs;

    for (int i=0; i<this->count(); ++i)
    {
        nrgs.append(this->operator()(i).energy(map));
    }

    return nrgs;
}

SireUnits::Dimension::MolarEnergy SelectorDihedral::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy SelectorDihedral::energy(const PropertyMap &map) const
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
    template class Mover<SelectorDihedral>;
}
