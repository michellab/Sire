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

using SireMM::detail::IDQuad;

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

SelectorImproper::SelectorImproper(const Improper &improper)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(improper)
{
    if (not improper.isEmpty())
    {
        auto atom0 = this->data().info().atomIdx(improper.ID().atom0());
        auto atom1 = this->data().info().atomIdx(improper.ID().atom1());
        auto atom2 = this->data().info().atomIdx(improper.ID().atom2());
        auto atom3 = this->data().info().atomIdx(improper.ID().atom3());


        if (atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or
            atom1 == atom2 or atom1 == atom3 or atom2 == atom3)
        {
            // cannot add Impropers to the same atom
            this->operator=(SelectorImproper());
        }
        else
        {
            this->imps.append(ImproperID(atom0, atom1, atom2, atom3));
        }
    }
}

inline
QSet<IDQuad> _to_int_set(const QList<ImproperID> &vals,
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
QList<ImproperID> _from_int_set(const QSet<IDQuad> &vals)
{
    QVector<IDQuad> v;
    v.reserve(vals.count());

    for (const auto &val : vals)
    {
        v.append(val);
    }

    std::sort(v.begin(), v.end());

    QList<ImproperID> l;
    l.reserve(v.count());

    for (const auto &val : v)
    {
        l.append(ImproperID(AtomIdx(val.atom0),
                            AtomIdx(val.atom1),
                            AtomIdx(val.atom2),
                            AtomIdx(val.atom3)));
    }

    return l;
}

QList<ImproperID> _get_impropers(const MoleculeData &moldata,
                                 const PropertyMap &map)
{
    QList<ImproperID> impropers;

    if (moldata.hasProperty(map["improper"]))
    {
        auto funcs = moldata.property(map["improper"]).asA<FourAtomFunctions>();

        bool have_connectivity = false;
        Connectivity c;

        if (moldata.hasProperty(map["connectivity"]))
        {
            try
            {
                c = moldata.property(map["connectivity"]).asA<Connectivity>();
                have_connectivity = true;
            }
            catch(...)
            {}
        }

        const auto &molinfo = moldata.info();

        for (const auto &func : funcs.potentials())
        {
            auto atomidx0 = molinfo.atomIdx(func.atom0());
            auto atomidx1 = molinfo.atomIdx(func.atom1());
            auto atomidx2 = molinfo.atomIdx(func.atom2());
            auto atomidx3 = molinfo.atomIdx(func.atom3());

            if (have_connectivity)
            {
                // we need to check the connectivity property too to find
                // the bond and angles so we can identify the central
                // atom in the improper
                if (c.areConnected(atomidx0, atomidx1) and
                    c.areConnected(atomidx2, atomidx1) and
                    c.areConnected(atomidx3, atomidx1))
                {
                    // this is the expected order :-)
                    impropers.append(ImproperID(atomidx0, atomidx1,
                                                atomidx2, atomidx3));
                }
                else if (c.areConnected(atomidx0, atomidx2) and
                         c.areConnected(atomidx1, atomidx2) and
                         c.areConnected(atomidx3, atomidx2))
                {
                    //this is the reverse order
                    impropers.append(ImproperID(atomidx0, atomidx2,
                                                atomidx1, atomidx3));
                }
                else if (c.areConnected(atomidx0, atomidx1) and
                         c.areConnected(atomidx0, atomidx2) and
                         c.areConnected(atomidx0, atomidx3))
                {
                    //the first atom is the central atom - strange
                    impropers.append(ImproperID(atomidx1, atomidx0,
                                                atomidx2, atomidx3));
                }
                else if (c.areConnected(atomidx3, atomidx0) and
                         c.areConnected(atomidx3, atomidx1) and
                         c.areConnected(atomidx3, atomidx2))
                {
                    //the last atom is the central atom - strange
                    impropers.append(ImproperID(atomidx0, atomidx3,
                                                atomidx1, atomidx2));
                }
                else
                {
                    auto i = ImproperID(atomidx0, atomidx1,
                                        atomidx2, atomidx3);

                    qDebug() << i.toString()
                             << "does not match the connectivity of the molecule.";

                    impropers.append(i);
                }
            }
            else
            {
                // we just have to trust this is right - normally from
                // amber files the third atom is the central atom
                impropers.append(ImproperID(atomidx0, atomidx2,
                                            atomidx1, atomidx3));
            }
        }
    }

    return impropers;
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    // do this so that the impropers are canonicalized and sorted
    auto impropers = _from_int_set(
                        _to_int_set(_get_impropers(this->data(), map),
                                    this->data().info()));

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

SelectorImproper::SelectorImproper(const MoleculeData &moldata,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>()
{
    this->operator=(SelectorImproper(Molecule(moldata), map));
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
QSet<IDQuad> _filter(const QSet<IDQuad> &impropers,
                     const QVector<quint32> &atoms,
                     int position)
{
    QSet<IDQuad> result;

    result.reserve(impropers.count());

    for (const auto &improper : impropers)
    {
        const auto atomidx = improper[position];

        for (const auto &atom : atoms)
        {
            if (atomidx == atom)
            {
                result.insert(improper);
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

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const ImproperID &improper,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(improper.atom0());
    auto atoms1 = mol.data().info().map(improper.atom1());
    auto atoms2 = mol.data().info().map(improper.atom2());
    auto atoms3 = mol.data().info().map(improper.atom3());

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

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                    this->data().info());

    auto impropers0123 = _filter(impropers, int_atoms0, 0);
    impropers0123 = _filter(impropers0123, int_atoms1, 1);
    impropers0123 = _filter(impropers0123, int_atoms2, 2);
    impropers0123 = _filter(impropers0123, int_atoms3, 3);

    imps = _from_int_set(impropers0123);
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const QList<ImproperID> &impropers,
                                   const SireBase::PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    if (impropers.count() == 1)
    {
        this->operator=(SelectorImproper(mol, impropers[0], map));
        return;
    }
    else if (not impropers.isEmpty())
    {
        auto all_impropers = _to_int_set(_get_impropers(this->data(), map),
                                         this->data().info());

        QSet<IDQuad> selected_impropers;
        selected_impropers.reserve(all_impropers.count());

        for (const auto &improper : impropers)
        {
            auto atoms0 = mol.data().info().map(improper.atom0());
            auto atoms1 = mol.data().info().map(improper.atom1());
            auto atoms2 = mol.data().info().map(improper.atom2());
            auto atoms3 = mol.data().info().map(improper.atom3());

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

            auto impropers0123 = _filter(all_impropers, int_atoms0, 0);
            impropers0123 = _filter(impropers0123, int_atoms1, 1);
            impropers0123 = _filter(impropers0123, int_atoms2, 2);
            impropers0123 = _filter(impropers0123, int_atoms3, 3);

            selected_impropers += impropers0123;

            if (selected_impropers.count() == all_impropers.count())
                break;
        }

        imps = _from_int_set(selected_impropers);
    }
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom, const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms = mol.data().info().map(atom);

    if (not mol.selectedAll())
    {
        const auto selection = mol.selection();
        atoms = _filter(atoms, selection);
    }

    auto int_atoms = _to_int(atoms);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers0 = _filter(impropers, int_atoms, 0);
    auto impropers1 = _filter(impropers, int_atoms, 1);
    auto impropers2 = _filter(impropers, int_atoms, 2);
    auto impropers3 = _filter(impropers, int_atoms, 3);

    imps = _from_int_set(impropers0 + impropers1 +
                         impropers2 + impropers3);
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);

    if (not mol.selectedAll())
    {
        const auto selection = mol.selection();
        atoms0 = _filter(atoms0, selection);
        atoms1 = _filter(atoms1, selection);
    }

    auto int_atoms0 = _to_int(atoms0);
    auto int_atoms1 = _to_int(atoms1);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers10 = _filter(impropers, int_atoms0, 1);
    impropers10 = _filter(impropers10, int_atoms1, 0);

    auto impropers12 = _filter(impropers, int_atoms0, 1);
    impropers12 = _filter(impropers12, int_atoms1, 2);

    auto impropers13 = _filter(impropers, int_atoms0, 1);
    impropers13 = _filter(impropers13, int_atoms1, 3);

    auto impropers01 = _filter(impropers, int_atoms0, 0);
    impropers01 = _filter(impropers01, int_atoms1, 1);

    auto impropers21 = _filter(impropers, int_atoms0, 2);
    impropers21 = _filter(impropers21, int_atoms1, 1);

    auto impropers31 = _filter(impropers, int_atoms0, 3);
    impropers31 = _filter(impropers31, int_atoms1, 1);

    imps = _from_int_set(impropers10 + impropers12 + impropers13 +
                         impropers01 + impropers21 + impropers31);
}

SelectorImproper::SelectorImproper(const MoleculeView &mol,
                                   const AtomID &atom0, const AtomID &atom1,
                                   const AtomID &atom2, const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(atom0);
    auto atoms1 = mol.data().info().map(atom1);
    auto atoms2 = mol.data().info().map(atom2);

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

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers012 = _filter(impropers, int_atoms0, 0);
    impropers012 = _filter(impropers012, int_atoms1, 1);
    impropers012 = _filter(impropers012, int_atoms2, 2);

    auto impropers013 = _filter(impropers, int_atoms0, 0);
    impropers013 = _filter(impropers013, int_atoms1, 1);
    impropers013 = _filter(impropers013, int_atoms2, 3);

    auto impropers213 = _filter(impropers, int_atoms0, 2);
    impropers213 = _filter(impropers213, int_atoms1, 1);
    impropers213 = _filter(impropers213, int_atoms2, 3);

    auto impropers210 = _filter(impropers, int_atoms0, 2);
    impropers210 = _filter(impropers210, int_atoms1, 1);
    impropers210 = _filter(impropers210, int_atoms2, 0);

    auto impropers310 = _filter(impropers, int_atoms0, 3);
    impropers310 = _filter(impropers310, int_atoms1, 1);
    impropers310 = _filter(impropers310, int_atoms2, 0);

    auto impropers312 = _filter(impropers, int_atoms0, 3);
    impropers312 = _filter(impropers312, int_atoms1, 1);
    impropers312 = _filter(impropers312, int_atoms2, 2);

    imps = _from_int_set(impropers012 + impropers013 + impropers213 +
                         impropers210 + impropers310 + impropers312);
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
    auto int_atoms = _to_int(atoms);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers0 = _filter(impropers, int_atoms, 0);
    auto impropers1 = _filter(impropers, int_atoms, 1);
    auto impropers2 = _filter(impropers, int_atoms, 2);
    auto impropers3 = _filter(impropers, int_atoms, 3);

    imps = _from_int_set(impropers0 + impropers1 +
                         impropers2 + impropers3);
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an improper from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    auto int_atoms0 = _to_int(atoms0);
    auto int_atoms1 = _to_int(atoms1);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers01 = _filter(impropers, int_atoms0, 0);
    impropers01 = _filter(impropers01, int_atoms1, 1);

    auto impropers21 = _filter(impropers, int_atoms0, 2);
    impropers21 = _filter(impropers21, int_atoms1, 1);

    auto impropers31 = _filter(impropers, int_atoms0, 3);
    impropers31 = _filter(impropers31, int_atoms1, 1);

    auto impropers10 = _filter(impropers, int_atoms0, 1);
    impropers10 = _filter(impropers10, int_atoms1, 0);

    auto impropers12 = _filter(impropers, int_atoms0, 1);
    impropers12 = _filter(impropers12, int_atoms1, 2);

    auto impropers13 = _filter(impropers, int_atoms0, 1);
    impropers13 = _filter(impropers13, int_atoms1, 3);

    imps = _from_int_set(impropers01 + impropers21 + impropers31 +
                         impropers10 + impropers12 + impropers13);
}

SelectorImproper::SelectorImproper(const Selector<Atom> &atoms0,
                                   const Selector<Atom> &atoms1,
                                   const Selector<Atom> &atoms2,
                                   const PropertyMap &map)
                 : ConcreteProperty<SelectorImproper, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and
             atoms0.isSameMolecule(atoms2)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an improper from atoms in the same molecule. "
            "%1, %2 and %3 are from different molecules (%4, %5 and %6)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString()), CODELOC);

    auto int_atoms0 = _to_int(atoms0);
    auto int_atoms1 = _to_int(atoms1);
    auto int_atoms2 = _to_int(atoms2);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers012 = _filter(impropers, int_atoms0, 0);
    impropers012 = _filter(impropers012, int_atoms1, 1);
    impropers012 = _filter(impropers012, int_atoms2, 2);

    auto impropers013 = _filter(impropers, int_atoms0, 0);
    impropers013 = _filter(impropers013, int_atoms1, 1);
    impropers013 = _filter(impropers013, int_atoms2, 3);

    auto impropers213 = _filter(impropers, int_atoms0, 2);
    impropers213 = _filter(impropers213, int_atoms1, 1);
    impropers213 = _filter(impropers213, int_atoms2, 3);

    auto impropers210 = _filter(impropers, int_atoms0, 2);
    impropers210 = _filter(impropers210, int_atoms1, 1);
    impropers210 = _filter(impropers210, int_atoms2, 0);

    auto impropers310 = _filter(impropers, int_atoms0, 3);
    impropers310 = _filter(impropers310, int_atoms1, 1);
    impropers310 = _filter(impropers310, int_atoms2, 0);

    auto impropers312 = _filter(impropers, int_atoms0, 3);
    impropers312 = _filter(impropers312, int_atoms1, 1);
    impropers312 = _filter(impropers312, int_atoms2, 2);

    imps = _from_int_set(impropers012 + impropers013 + impropers213 +
                         impropers210 + impropers310 + impropers312);
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
            "You can only create an improper from atoms in the same molecule. "
            "%1, %2, %3 and %4 are from different molecules (%5, %6, %7 and %8)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString()).arg(atoms3.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString())
                .arg(atoms3.molecule().toString()), CODELOC);

    auto int_atoms0 = _to_int(atoms0);
    auto int_atoms1 = _to_int(atoms1);
    auto int_atoms2 = _to_int(atoms2);
    auto int_atoms3 = _to_int(atoms3);

    auto impropers = _to_int_set(_get_impropers(this->data(), map),
                                 this->data().info());

    auto impropers0123 = _filter(impropers, int_atoms0, 0);
    impropers0123 = _filter(impropers0123, int_atoms1, 1);
    impropers0123 = _filter(impropers0123, int_atoms2, 2);
    impropers0123 = _filter(impropers0123, int_atoms3, 3);

    imps = _from_int_set(impropers0123);
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

QList<SireUnits::Dimension::Angle> SelectorImproper::measures(const PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireUnits::Dimension::Angle> SelectorImproper::measures() const
{
    return this->sizes();
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
