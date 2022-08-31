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

using SireMM::detail::IDQuad;

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

SelectorDihedral::SelectorDihedral(const Dihedral &dihedral)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(dihedral)
{
    if (not dihedral.isEmpty())
    {
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
        {
            // cannot add Dihedrals to the same atom
            this->operator=(SelectorDihedral());
        }
        else
        {
            this->dihs.append(DihedralID(atom0, atom1, atom2, atom3));
        }
    }
}

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

SelectorDihedral::SelectorDihedral(const MoleculeData &moldata,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>()
{
    this->operator=(SelectorDihedral(Molecule(moldata), map));
}

inline
QSet<IDQuad> _to_int_set(const QList<DihedralID> &vals,
                         const MoleculeInfoData &molinfo)
{
    QSet<IDQuad> s;
    s.reserve(vals.count());

    for (const auto &val : vals)
    {
        s.insert(IDQuad(molinfo.atomIdx(val[0]).value(),
                        molinfo.atomIdx(val[1]).value(),
                        molinfo.atomIdx(val[2]).value(),
                        molinfo.atomIdx(val[3]).value()));
    }

    return s;
}

inline
QList<DihedralID> _from_int_set(const QSet<IDQuad> &vals)
{
    QVector<IDQuad> v;
    v.reserve(vals.count());

    for (const auto &val : vals)
    {
        v.append(val);
    }

    std::sort(v.begin(), v.end());

    QList<DihedralID> l;
    l.reserve(v.count());

    for (const auto &val : v)
    {
        l.append(DihedralID(AtomIdx(val.atom0),
                            AtomIdx(val.atom1),
                            AtomIdx(val.atom2),
                            AtomIdx(val.atom3)));
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
QSet<IDQuad> _filter(const QSet<IDQuad> &dihedrals,
                     const QVector<quint32> &atoms,
                     int position)
{
    QSet<IDQuad> result;

    result.reserve(dihedrals.count());

    for (const auto &dihedral : dihedrals)
    {
        const auto atomidx = dihedral[position];

        for (const auto &atom : atoms)
        {
            if (atomidx == atom)
            {
                result.insert(dihedral);
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
        if (not mol.selectedAll())
        {
            const auto selection = mol.selection();
            atoms0 = _filter(atoms0, selection);
            atoms1 = _filter(atoms1, selection);
            atoms2 = _filter(atoms2, selection);
            atoms3 = _filter(atoms3, selection);
        }

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);
        auto int_atoms2 = _to_int(atoms2);
        auto int_atoms3 = _to_int(atoms3);

        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals0123 = _filter(dihedrals, int_atoms0, 0);
        dihedrals0123 = _filter(dihedrals0123, int_atoms1, 1);
        dihedrals0123 = _filter(dihedrals0123, int_atoms2, 2);
        dihedrals0123 = _filter(dihedrals0123, int_atoms3, 3);

        auto dihedrals3210 = _filter(dihedrals, int_atoms0, 3);
        dihedrals3210 = _filter(dihedrals3210, int_atoms1, 2);
        dihedrals3210 = _filter(dihedrals3210, int_atoms2, 1);
        dihedrals3210 = _filter(dihedrals3210, int_atoms3, 0);

        dihs = _from_int_set(dihedrals0123 + dihedrals3210);
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const QList<DihedralID> &dihedrals,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    if (dihedrals.count() == 1)
    {
        this->operator=(SelectorDihedral(mol, dihedrals[0], map));
        return;
    }
    else if (mol.data().hasProperty(map["connectivity"]) and not dihedrals.isEmpty())
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto all_dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        QSet<IDQuad> selected_dihedrals;
        selected_dihedrals.reserve(all_dihedrals.count());

        for (const auto &dihedral : dihedrals)
        {
            auto atoms0 = mol.data().info().map(dihedral.atom0());
            auto atoms1 = mol.data().info().map(dihedral.atom1());
            auto atoms2 = mol.data().info().map(dihedral.atom2());
            auto atoms3 = mol.data().info().map(dihedral.atom3());

            if (not mol.selectedAll())
            {
                const auto selection = mol.selection();
                atoms0 = _filter(atoms0, selection);
                atoms1 = _filter(atoms1, selection);
                atoms2 = _filter(atoms2, selection);
                atoms3 = _filter(atoms3, selection);
            }

            auto int_atoms0 = _to_int(atoms0);
            auto int_atoms1 = _to_int(atoms1);
            auto int_atoms2 = _to_int(atoms2);
            auto int_atoms3 = _to_int(atoms3);

            auto dihedrals0123 = _filter(all_dihedrals, int_atoms0, 0);
            dihedrals0123 = _filter(dihedrals0123, int_atoms1, 1);
            dihedrals0123 = _filter(dihedrals0123, int_atoms2, 2);
            dihedrals0123 = _filter(dihedrals0123, int_atoms3, 3);

            auto dihedrals3210 = _filter(all_dihedrals, int_atoms0, 3);
            dihedrals3210 = _filter(dihedrals3210, int_atoms1, 2);
            dihedrals3210 = _filter(dihedrals3210, int_atoms2, 1);
            dihedrals3210 = _filter(dihedrals3210, int_atoms3, 0);

            selected_dihedrals += dihedrals0123;
            selected_dihedrals += dihedrals3210;

            if (selected_dihedrals.count() == all_dihedrals.count())
                break;
        }

        dihs = _from_int_set(selected_dihedrals);
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
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

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals0 = _filter(dihedrals, int_atoms, 0);
        auto dihedrals1 = _filter(dihedrals, int_atoms, 1);
        auto dihedrals2 = _filter(dihedrals, int_atoms, 2);
        auto dihedrals3 = _filter(dihedrals, int_atoms, 3);

        dihs = _from_int_set(dihedrals0 + dihedrals1 +
                             dihedrals2 + dihedrals3);
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);

    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        if (not mol.selectedAll())
        {
            const auto selection = mol.selection();
            atoms0 = _filter(atoms0, selection);
            atoms1 = _filter(atoms1, selection);
        }

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals01 = _filter(dihedrals, int_atoms0, 0);
        dihedrals01 = _filter(dihedrals01, int_atoms1, 1);

        auto dihedrals12 = _filter(dihedrals, int_atoms0, 1);
        dihedrals12 = _filter(dihedrals12, int_atoms1, 2);

        auto dihedrals23 = _filter(dihedrals, int_atoms0, 2);
        dihedrals23 = _filter(dihedrals23, int_atoms1, 3);

        auto dihedrals32 = _filter(dihedrals, int_atoms0, 3);
        dihedrals32 = _filter(dihedrals32, int_atoms1, 2);

        auto dihedrals21 = _filter(dihedrals, int_atoms0, 2);
        dihedrals21 = _filter(dihedrals21, int_atoms1, 1);

        auto dihedrals10 = _filter(dihedrals, int_atoms0, 1);
        dihedrals10 = _filter(dihedrals, int_atoms1, 0);

        dihs = _from_int_set(dihedrals01 + dihedrals12 + dihedrals23 +
                             dihedrals32 + dihedrals21 + dihedrals10);
    }
}

SelectorDihedral::SelectorDihedral(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);
    auto atoms2 = mol.data().info().map(atom2);

    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        if (not mol.selectedAll())
        {
            const auto selection = mol.selection();
            atoms0 = _filter(atoms0, selection);
            atoms1 = _filter(atoms1, selection);
            atoms2 = _filter(atoms2, selection);
        }

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);
        auto int_atoms2 = _to_int(atoms2);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals012 = _filter(dihedrals, int_atoms0, 0);
        dihedrals012 = _filter(dihedrals012, int_atoms1, 1);
        dihedrals012 = _filter(dihedrals012, int_atoms2, 2);

        auto dihedrals123 = _filter(dihedrals, int_atoms0, 1);
        dihedrals123 = _filter(dihedrals123, int_atoms1, 2);
        dihedrals123 = _filter(dihedrals123, int_atoms2, 3);

        auto dihedrals321 = _filter(dihedrals, int_atoms0, 3);
        dihedrals321 = _filter(dihedrals321, int_atoms1, 2);
        dihedrals321 = _filter(dihedrals321, int_atoms2, 1);

        auto dihedrals210 = _filter(dihedrals, int_atoms0, 2);
        dihedrals210 = _filter(dihedrals210, int_atoms1, 1);
        dihedrals210 = _filter(dihedrals210, int_atoms2, 0);

        dihs = _from_int_set(dihedrals012 + dihedrals123 +
                             dihedrals321 + dihedrals210);
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
    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms = _to_int(atoms);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals0 = _filter(dihedrals, int_atoms, 0);
        auto dihedrals1 = _filter(dihedrals, int_atoms, 1);
        auto dihedrals2 = _filter(dihedrals, int_atoms, 2);
        auto dihedrals3 = _filter(dihedrals, int_atoms, 3);

        dihs = _from_int_set(dihedrals0 + dihedrals1 +
                             dihedrals2 + dihedrals3);
    }
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a dihedral from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals01 = _filter(dihedrals, int_atoms0, 0);
        dihedrals01 = _filter(dihedrals01, int_atoms1, 1);

        auto dihedrals12 = _filter(dihedrals, int_atoms0, 1);
        dihedrals12 = _filter(dihedrals12, int_atoms1, 2);

        auto dihedrals23 = _filter(dihedrals, int_atoms0, 2);
        dihedrals23 = _filter(dihedrals23, int_atoms1, 3);

        auto dihedrals32 = _filter(dihedrals, int_atoms0, 3);
        dihedrals32 = _filter(dihedrals32, int_atoms1, 2);

        auto dihedrals21 = _filter(dihedrals, int_atoms0, 2);
        dihedrals21 = _filter(dihedrals21, int_atoms1, 1);

        auto dihedrals10 = _filter(dihedrals, int_atoms0, 1);
        dihedrals10 = _filter(dihedrals10, int_atoms1, 0);

        dihs = _from_int_set(dihedrals01 + dihedrals12 + dihedrals23 +
                             dihedrals32 + dihedrals21 + dihedrals10);
    }
}

SelectorDihedral::SelectorDihedral(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorDihedral, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and
             atoms0.isSameMolecule(atoms2)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a dihedral from atoms in the same molecule. "
            "%1, %2 and %3 are from different molecules (%4, %5 and %6)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString()), CODELOC);

    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);
        auto int_atoms2 = _to_int(atoms2);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals012 = _filter(dihedrals, int_atoms0, 0);
        dihedrals012 = _filter(dihedrals012, int_atoms1, 1);
        dihedrals012 = _filter(dihedrals012, int_atoms2, 2);

        auto dihedrals123 = _filter(dihedrals, int_atoms0, 1);
        dihedrals123 = _filter(dihedrals123, int_atoms1, 2);
        dihedrals123 = _filter(dihedrals123, int_atoms2, 3);

        auto dihedrals321 = _filter(dihedrals, int_atoms0, 3);
        dihedrals321 = _filter(dihedrals321, int_atoms1, 2);
        dihedrals321 = _filter(dihedrals321, int_atoms2, 1);

        auto dihedrals210 = _filter(dihedrals, int_atoms0, 2);
        dihedrals210 = _filter(dihedrals210, int_atoms1, 1);
        dihedrals210 = _filter(dihedrals210, int_atoms2, 0);

        dihs = _from_int_set(dihedrals012 + dihedrals123 +
                             dihedrals321 + dihedrals210);
    }
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
            "You can only create a dihedral from atoms in the same molecule. "
            "%1, %2, %3 and %4 are from different molecules (%5, %6, %7 and %8)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString()).arg(atoms3.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString())
                .arg(atoms3.molecule().toString()), CODELOC);

    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);
        auto int_atoms2 = _to_int(atoms2);
        auto int_atoms3 = _to_int(atoms3);

        auto dihedrals = _to_int_set(c.getDihedrals(), this->data().info());

        auto dihedrals0123 = _filter(dihedrals, int_atoms0, 0);
        dihedrals0123 = _filter(dihedrals0123, int_atoms1, 1);
        dihedrals0123 = _filter(dihedrals0123, int_atoms2, 2);
        dihedrals0123 = _filter(dihedrals0123, int_atoms3, 3);

        auto dihedrals3210 = _filter(dihedrals, int_atoms0, 3);
        dihedrals3210 = _filter(dihedrals3210, int_atoms1, 2);
        dihedrals3210 = _filter(dihedrals3210, int_atoms2, 1);
        dihedrals3210 = _filter(dihedrals3210, int_atoms3, 0);

        dihs = _from_int_set(dihedrals0123 + dihedrals3210);
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

    auto ret_dihs = _to_int_set(this->dihs, this->data().info());
    ret_dihs.intersect(_to_int_set(other.dihs, other.data().info()));

    SelectorDihedral ret(*this);
    ret.dihs = _from_int_set(ret_dihs);

    return ret;
}

SelectorDihedral SelectorDihedral::invert(const PropertyMap &map) const
{
    auto s = SelectorDihedral(this->molecule(), map);

    auto ret_dihs = _to_int_set(s.dihs, s.data().info());

    ret_dihs.subtract(_to_int_set(this->dihs, this->data().info()));

    SelectorDihedral ret(*this);
    ret.dihs = _from_int_set(ret_dihs);

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
