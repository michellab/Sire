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

#include "selectorangle.h"

#include "threeatomfunctions.h"

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

static const RegisterMetaType<SelectorAngle> r_sangle;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const SelectorAngle &sangle)
{
    writeHeader(ds, r_sangle, 1);
    SharedDataStream sds(ds);
    sds << sangle.angs << static_cast<const MoleculeView&>(sangle);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, SelectorAngle &sangle)
{
    auto v = readHeader(ds, r_sangle);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> sangle.angs >> static_cast<MoleculeView&>(sangle);
    }
    else
        throw SireStream::version_error(v, "1", r_sangle, CODELOC);

    return ds;
}

SelectorAngle::SelectorAngle() : ConcreteProperty<SelectorAngle, MoleculeView>()
{}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto angles = c.getAngles();

        if (mol.selectedAll())
        {
            angs = angles;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &angle : angles)
            {
                if (s.selected(angle.atom0()) and
                    s.selected(angle.atom1()) and
                    s.selected(angle.atom2()))
                {
                    angs.append(angle);
                }
            }
        }
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const QList<AngleID> &angles)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    const auto s = mol.selection();

    for (const auto &angle : angles)
    {
        AngleID a(mol.data().info().atomIdx(angle.atom0()),
                  mol.data().info().atomIdx(angle.atom1()),
                  mol.data().info().atomIdx(angle.atom2()));

        if (s.selected(a.atom0()) and
            s.selected(a.atom1()) and
            s.selected(a.atom2()))
        {
            angs.append(a);
        }
    }
}

SelectorAngle::SelectorAngle(const MoleculeData &moldata,
                             const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(Molecule(moldata), map));
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AngleID &angle, const PropertyMap &map)
             : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(angle.atom0());
    auto atoms1 = mol.data().info().map(angle.atom1());
    auto atoms2 = mol.data().info().map(angle.atom2());

    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<AngleID> seen_angles;

        QList<AngleID> angles;

        for (const auto &atom0 : atoms0)
        {
            for (const auto &atom1 : atoms1)
            {
                for (const auto &atom2 : atoms2)
                {
                    auto atomidx0 = atom0;
                    auto atomidx1 = atom1;
                    auto atomidx2 = atom2;

                    if (atomidx0 > atomidx2)
                    {
                        qSwap(atomidx0, atomidx2);
                    }

                    AngleID a(atomidx0, atomidx1, atomidx2);

                    if (atomidx0 != atomidx1 and
                        atomidx0 != atomidx2 and
                        atomidx1 != atomidx2 and
                        c.areConnected(atomidx0, atomidx1) and
                        c.areConnected(atomidx1, atomidx2) and
                        not seen_angles.contains(a))
                    {
                        seen_angles.insert(a);
                        angles.append(a);
                    }
                }
            }
        }

        if (mol.selectedAll())
        {
            angs = angles;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &angle : angles)
            {
                if (s.selected(angle.atom0()) and
                    s.selected(angle.atom1()) and
                    s.selected(angle.atom2()))
                {
                    angs.append(angle);
                }
            }
        }
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AtomID &atom, const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto angles = c.getAngles(atom);

        if (mol.selectedAll())
        {
            angs = angles;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &angle : angles)
            {
                if (s.selected(angle.atom0()) and
                    s.selected(angle.atom1()) and
                    s.selected(angle.atom2()))
                {
                    angs.append(angle);
                }
            }
        }
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AtomID &atom0, const AtomID &atom1,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto angles = c.getAngles(atom0, atom1);

        if (mol.selectedAll())
        {
            angs = angles;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &angle : angles)
            {
                if (s.selected(angle.atom0()) and
                    s.selected(angle.atom1()) and
                    s.selected(angle.atom2()))
                {
                    angs.append(angle);
                }
            }
        }
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AtomID &atom0, const AtomID &atom1,
                             const AtomID &atom2,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(mol, AngleID(atom0, atom1, atom2), map));
}

SelectorAngle::SelectorAngle(const MoleculeData &mol,
                             const AtomID &atom, const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(Molecule(mol), atom, map));
}

SelectorAngle::SelectorAngle(const MoleculeData &mol,
                             const AtomID &atom0, const AtomID &atom1,
                             const PropertyMap &map)
             : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(Molecule(mol), atom0, atom1, map));
}

SelectorAngle::SelectorAngle(const MoleculeData &mol,
                             const AtomID &atom0, const AtomID &atom1,
                             const AtomID &atom2,
                             const PropertyMap &map)
             : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(Molecule(mol), atom0, atom1, atom2, map));
}

SelectorAngle::SelectorAngle(const Selector<Atom> &atoms,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(atoms)
{
    if (atoms.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<AngleID> seen_angs;

        QList<AngleID> angles;

        for (int i=0; i<atoms.count(); ++i)
        {
            for (const auto &a : c.getAngles(atoms(i).index()))
            {
                auto atomidx0 = atoms.data().info().atomIdx(a.atom0());
                auto atomidx1 = atoms.data().info().atomIdx(a.atom1());
                auto atomidx2 = atoms.data().info().atomIdx(a.atom2());

                if (atomidx0 > atomidx2)
                {
                    qSwap(atomidx0, atomidx2);
                }

                AngleID ang(atomidx0, atomidx1, atomidx2);

                if (atomidx0 != atomidx1 and
                    atomidx0 != atomidx2 and
                    atomidx1 != atomidx2 and
                    not seen_angs.contains(ang))
                {
                    seen_angs.insert(ang);
                    angles.append(ang);
                }
            }
        }

        angs = angles;
    }
}

SelectorAngle::SelectorAngle(const Selector<Atom> &atoms0,
                             const Selector<Atom> &atoms1,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an Angle from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        QSet<AngleID> seen_angs;

        QList<AngleID> angles;

        for (int i=0; i<atoms0.count(); ++i)
        {
            for (int j=0; j<atoms1.count(); ++j)
            {
                auto atomidx0 = atoms0(i).index();
                auto atomidx1 = atoms1(j).index();

                auto angs = c.getAngles(atomidx0, atomidx1);

                for (const auto &ang : angs)
                {
                    if (not seen_angs.contains(ang))
                    {
                        seen_angs.insert(ang);
                        angles.append(ang);
                    }
                }
            }
        }

        angs = angles;
    }
}

bool _contains(const AngleID &ang, const AtomIdx &atom)
{
    return ang.atom0() == atom or
           ang.atom1() == atom or
           ang.atom2() == atom;
}

SelectorAngle::SelectorAngle(const Selector<Atom> &atoms0,
                             const Selector<Atom> &atoms1,
                             const Selector<Atom> &atoms2,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(atoms0)
{
    if (not (atoms0.isSameMolecule(atoms1) and
             atoms0.isSameMolecule(atoms2)))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an Angle from atoms in the same molecule. "
            "%1, %2 and %3 are from different molecules (%4, %5 and %6)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms2.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString())
                .arg(atoms2.molecule().toString()), CODELOC);

    if (atoms0.data().hasProperty(map["connectivity"]))
    {
        auto c = atoms0.data().property(map["connectivity"]).asA<Connectivity>();

        bool found = false;

        for (const auto &angle : c.getAngles())
        {
            for (int i=0; i<atoms0.count(); ++i)
            {
                auto atomidx0 = atoms0(i).index();

                if (not _contains(angle, atomidx0))
                    break;

                for (int j=0; j<atoms1.count(); ++j)
                {
                    auto atomidx1 = atoms1(j).index();

                    if (atomidx1 == atomidx0 or not _contains(angle, atomidx1))
                        break;

                    for (int k=0; k<atoms2.count(); ++k)
                    {
                        auto atomidx2 = atoms2(k).index();

                        if (atomidx2 == atomidx1 or atomidx2 == atomidx0 or
                            not _contains(angle, atomidx2))
                            break;

                        if (atomidx0 > atomidx2)
                        {
                            qSwap(atomidx0, atomidx2);
                        }

                        angs.append(AngleID(atomidx0, atomidx1, atomidx2));
                        found = true;
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

SelectorAngle::SelectorAngle(const SelectorAngle &other)
              : ConcreteProperty<SelectorAngle, MoleculeView>(other),
                angs(other.angs)
{}

SelectorAngle::~SelectorAngle()
{}

const char* SelectorAngle::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorAngle>());
}

SelectorAngle& SelectorAngle::operator=(const SelectorAngle &other)
{
    if (this != &other)
    {
        angs = other.angs;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorAngle::operator==(const SelectorAngle &other) const
{
    return angs == other.angs and MoleculeView::operator==(other);
}

bool SelectorAngle::operator!=(const SelectorAngle &other) const
{
    return not operator==(other);
}

int SelectorAngle::count() const
{
    return this->angs.count();
}

int SelectorAngle::size() const
{
    return this->count();
}

int SelectorAngle::nViews() const
{
    return this->count();
}

QString SelectorAngle::toString() const
{
    if (this->isNull())
        return QObject::tr("SelectorAngle::null");
    else if (this->isEmpty())
        return QObject::tr("SelectorAngle::empty");

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

    return QObject::tr("SelectorAngle( size=%1\n%2\n)")
                .arg(this->count()).arg(parts.join("\n"));
}

SelectorAngle SelectorAngle::add(const Angle &angle) const
{
    if (angle.isNull())
        return *this;

    if (this->isEmpty())
    {
        return SelectorAngle(angle);
    }

    if (angle.data().number() != this->data().number())
    {
        throw SireError::incompatible_error(QObject::tr(
            "You cannot add Angles from a different molecule (%1) to "
            "a set of Angles from molecule %2.")
                .arg(angle.data().number())
                .arg(this->data().number()),
                    CODELOC);
    }

    auto atom0 = this->data().info().atomIdx(angle.ID().atom0());
    auto atom1 = this->data().info().atomIdx(angle.ID().atom1());
    auto atom2 = this->data().info().atomIdx(angle.ID().atom2());

    if (atom0 > atom2)
        qSwap(atom0, atom2);

    if (atom0 == atom1 or atom0 == atom2 or atom1 == atom2)
        // cannot add Angles to the same atom
        return *this;

    SelectorAngle ret(*this);

    ret.angs.append(AngleID(atom0, atom1, atom2));

    return ret;
}

MolViewPtr SelectorAngle::operator[](int i) const
{
    return this->operator()(i);
}

MolViewPtr SelectorAngle::operator[](const SireBase::Slice &slice) const
{
    return this->operator()(slice);
}

MolViewPtr SelectorAngle::operator[](const QList<qint64> &idxs) const
{
    return this->operator()(idxs);
}

MolViewPtr SelectorAngle::operator[](const AngleID &angle) const
{
    return this->operator()(angle);
}

SireMM::Angle SelectorAngle::operator()(int i) const
{
    auto ang = angs.at(Index(i).map(angs.count()));
    return SireMM::Angle(this->data(), ang);
}

SelectorAngle SelectorAngle::operator()(const SireBase::Slice &slice) const
{
    SelectorAngle ret(*this);
    ret.angs.clear();

    for (auto it = slice.begin(angs.count()); not it.atEnd(); it.next())
    {
        ret.angs.append(this->angs.at(it.value()));
    }

    return ret;
}

SelectorAngle SelectorAngle::operator()(const QList<qint64> &idxs) const
{
    SelectorAngle ret(*this);
    ret.angs.clear();

    for (const auto &idx : idxs)
    {
        ret.angs.append(this->angs.at(Index(idx).map(this->angs.count())));
    }

    return ret;
}

SelectorAngle SelectorAngle::operator()(int i, int j) const
{
    i = Index(i).map(angs.count());
    j = Index(j).map(angs.count());

    SelectorAngle ret(*this);
    ret.angs.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret.angs.append(this->angs.at(i));
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret.angs.append(this->angs.at(i));
        }
    }

    return ret;
}

SelectorAngle SelectorAngle::operator()(const AngleID &angle) const
{
    auto atom0s = this->data().info().map(angle.atom0());
    auto atom1s = this->data().info().map(angle.atom1());
    auto atom2s = this->data().info().map(angle.atom2());

    SelectorAngle ret(*this);
    ret.angs.clear();

    for (const auto &atom0 : atom0s)
    {
        for (const auto &atom1 : atom1s)
        {
            for (const auto &atom2 : atom2s)
            {
                auto a0 = atom0;
                auto a1 = atom1;
                auto a2 = atom2;

                if (a0 > a2)
                    qSwap(a0, a2);

                AngleID angle(a0, a1, a2);

                for (const auto &a : angs)
                {
                    if (a == angle)
                        ret.angs.append(a);
                }
            }
        }
    }

    return ret;
}

MolViewPtr SelectorAngle::toSelector() const
{
    return MolViewPtr(*this);
}

QList<MolViewPtr> SelectorAngle::toList() const
{
    QList<MolViewPtr> l;
    l.reserve(angs.count());

    auto d = this->data();

    for (const auto &ang : angs)
    {
        l.append(MolViewPtr(new Angle(d, ang)));
    }

    return l;
}

SelectorAngle SelectorAngle::add(const SelectorAngle &other) const
{
    if (this->isEmpty())
        return other;
    else if (other.isEmpty())
        return *this;

    MoleculeView::assertSameMolecule(other);

    SelectorAngle ret(*this);

    for (const auto &angle : other.angs)
    {
        if (not this->angs.contains(angle))
        {
            ret.angs.append(angle);
        }
    }

    return ret;
}

SelectorAngle SelectorAngle::intersection(const SelectorAngle &other) const
{
    if (this->isEmpty() or other.isEmpty())
    {
        return SelectorAngle();
    }

    MoleculeView::assertSameMolecule(other);

    SelectorAngle ret(*this);
    ret.angs.clear();

    for (const auto &angle : this->angs)
    {
        if (ret.angs.contains(angle))
        {
            ret.angs.append(angle);
        }
    }

    return ret;
}

SelectorAngle SelectorAngle::invert(const PropertyMap &map) const
{
    auto s = SelectorAngle(this->molecule(), map);

    SelectorAngle ret(*this);
    ret.angs.clear();

    for (const auto &angle : s.angs)
    {
        if (not this->angs.contains(angle))
        {
            ret.angs.append(angle);
        }
    }

    return ret;
}

SelectorAngle SelectorAngle::invert() const
{
    return this->invert(PropertyMap());
}

QList<AngleID> SelectorAngle::IDs() const
{
    return angs;
}

bool SelectorAngle::isEmpty() const
{
    return this->angs.isEmpty();
}

bool SelectorAngle::selectedAll() const
{
    return this->selection().selectedAll();
}

AtomSelection SelectorAngle::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    for (const auto &ang : angs)
    {
        s = s.select(ang.atom0());
        s = s.select(ang.atom1());
        s = s.select(ang.atom2());
    }

    return s;
}

bool SelectorAngle::hasProperty(const PropertyName &key) const
{
    for (int i=0; i<this->count(); ++i)
    {
        if (this->operator()(i).hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorAngle::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool SelectorAngle::hasMetadata(const PropertyName &key,
                                const PropertyName &metakey) const
{
    return false;
}

QStringList SelectorAngle::propertyKeys() const
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

QStringList SelectorAngle::metadataKeys() const
{
    return QStringList();
}

QStringList SelectorAngle::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

QList<Properties> SelectorAngle::properties() const
{
    QList<Properties> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).properties());
    }

    return props;
}

Mover<SelectorAngle> SelectorAngle::move() const
{
    return Mover<SelectorAngle>(*this);
}

Evaluator SelectorAngle::evaluate() const
{
    return Evaluator(*this);
}

QList<PropertyPtr> SelectorAngle::property(const PropertyName &key) const
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
            "None of the Angles in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<PropertyPtr> SelectorAngle::property(const PropertyName &key,
                                           const Property &default_value) const
{
    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).property(key, default_value));
    }

    return props;
}

QList<SireUnits::Dimension::Angle> SelectorAngle::sizes(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::Angle> s;

    for (int i=0; i<this->count(); ++i)
    {
        s.append(this->operator()(i).size(map));
    }

    return s;
}

QList<SireUnits::Dimension::Angle> SelectorAngle::sizes() const
{
    return this->sizes(PropertyMap());
}

QList<Expression> SelectorAngle::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<Expression> SelectorAngle::potentials(const PropertyMap &map) const
{
    QList<Expression> p;

    for (int i=0; i<this->count(); ++i)
    {
        p.append(this->operator()(i).potential(map));
    }

    return p;
}

QList<SireUnits::Dimension::MolarEnergy> SelectorAngle::energies() const
{
    return this->energies(PropertyMap());
}

QList<SireUnits::Dimension::MolarEnergy> SelectorAngle::energies(const PropertyMap &map) const
{
    QList<SireUnits::Dimension::MolarEnergy> nrgs;

    for (int i=0; i<this->count(); ++i)
    {
        nrgs.append(this->operator()(i).energy(map));
    }

    return nrgs;
}

SireUnits::Dimension::MolarEnergy SelectorAngle::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy SelectorAngle::energy(const PropertyMap &map) const
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
    template class Mover<SelectorAngle>;
}
