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

#include "SireMM/threeatomfunctions.h"

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

using SireMM::detail::IDTriple;

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

SelectorAngle::SelectorAngle(const MoleculeData &moldata,
                             const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorAngle, MoleculeView>()
{
    this->operator=(SelectorAngle(Molecule(moldata), map));
}

inline
QSet<IDTriple> _to_int_set(const QList<AngleID> &vals,
                           const MoleculeInfoData &molinfo)
{
    QSet<IDTriple> s;
    s.reserve(vals.count());

    for (const auto &val : vals)
    {
        s.insert(IDTriple(molinfo.atomIdx(val[0]).value(),
                          molinfo.atomIdx(val[1]).value(),
                          molinfo.atomIdx(val[2]).value()));
    }

    return s;
}

inline
QList<AngleID> _from_int_set(const QSet<IDTriple> &vals)
{
    QVector<IDTriple> v;
    v.reserve(vals.count());

    for (const auto &val : vals)
    {
        v.append(val);
    }

    std::sort(v.begin(), v.end());

    QList<AngleID> l;
    l.reserve(v.count());

    for (const auto &val : v)
    {
        l.append(AngleID(AtomIdx(val.atom0),
                         AtomIdx(val.atom1),
                         AtomIdx(val.atom2)));
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
QSet<IDTriple> _filter(const QSet<IDTriple> &angles,
                       const QVector<quint32> &atoms,
                       int position)
{
    QSet<IDTriple> result;

    result.reserve(angles.count());

    for (const auto &angle : angles)
    {
        const auto atomidx = angle[position];

        for (const auto &atom : atoms)
        {
            if (atomidx == atom)
            {
                result.insert(angle);
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

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AngleID &angle, const PropertyMap &map)
             : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    auto atoms0 = mol.data().info().map(angle.atom0());
    auto atoms1 = mol.data().info().map(angle.atom1());
    auto atoms2 = mol.data().info().map(angle.atom2());

    if (mol.data().hasProperty(map["connectivity"]))
    {
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

        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto left_angles = _filter(angles, int_atoms0, 0);
        left_angles = _filter(left_angles, int_atoms1, 1);
        left_angles = _filter(left_angles, int_atoms2, 2);

        auto right_angles = _filter(angles, int_atoms0, 2);
        right_angles = _filter(right_angles, int_atoms1, 1);
        right_angles = _filter(right_angles, int_atoms2, 0);

        angs = _from_int_set(left_angles + right_angles);
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const QList<AngleID> &angles,
                             const SireBase::PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
{
    if (angles.count() == 1)
    {
        this->operator=(SelectorAngle(mol, angles[0], map));
        return;
    }
    else if (mol.data().hasProperty(map["connectivity"]) and not angles.isEmpty())
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto all_angles = _to_int_set(c.getAngles(), this->data().info());

        QSet<IDTriple> selected_angles;
        selected_angles.reserve(all_angles.count());

        for (const auto &angle : angles)
        {
            auto atoms0 = mol.data().info().map(angle.atom0());
            auto atoms1 = mol.data().info().map(angle.atom1());
            auto atoms2 = mol.data().info().map(angle.atom2());

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

            auto left_angles = _filter(all_angles, int_atoms0, 0);
            left_angles = _filter(left_angles, int_atoms1, 1);
            left_angles = _filter(left_angles, int_atoms2, 2);

            auto right_angles = _filter(all_angles, int_atoms0, 2);
            right_angles = _filter(right_angles, int_atoms1, 1);
            right_angles = _filter(right_angles, int_atoms2, 0);

            selected_angles += left_angles;
            selected_angles += right_angles;

            if (selected_angles.count() == all_angles.count())
                break;
        }

        angs = _from_int_set(selected_angles);
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AtomID &atom, const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
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

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto angles0 = _filter(angles, int_atoms, 0);
        auto angles1 = _filter(angles, int_atoms, 1);
        auto angles2 = _filter(angles, int_atoms, 2);

        angs = _from_int_set(angles0 + angles1 + angles2);
    }
}

SelectorAngle::SelectorAngle(const MoleculeView &mol,
                             const AtomID &atom0, const AtomID &atom1,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(mol)
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

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto angles01 = _filter(angles, int_atoms0, 0);
        angles01 = _filter(angles01, int_atoms1, 1);

        auto angles12 = _filter(angles, int_atoms0, 1);
        angles12 = _filter(angles12, int_atoms1, 2);

        auto angles21 = _filter(angles, int_atoms0, 2);
        angles21 = _filter(angles21, int_atoms1, 1);

        auto angles10 = _filter(angles, int_atoms0, 1);
        angles10 = _filter(angles, int_atoms1, 0);

        angs = _from_int_set(angles01 + angles12 + angles21 + angles10);
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
    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms = _to_int(atoms);

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto angles0 = _filter(angles, int_atoms, 0);
        auto angles1 = _filter(angles, int_atoms, 1);
        auto angles2 = _filter(angles, int_atoms, 2);

        angs = _from_int_set(angles0 + angles1 + angles2);
    }
}

SelectorAngle::SelectorAngle(const Selector<Atom> &atoms0,
                             const Selector<Atom> &atoms1,
                             const PropertyMap &map)
              : ConcreteProperty<SelectorAngle, MoleculeView>(atoms0)
{
    if (not atoms0.isSameMolecule(atoms1))
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an angle from atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atoms0.toString()).arg(atoms1.toString())
                .arg(atoms0.molecule().toString())
                .arg(atoms1.molecule().toString()), CODELOC);

    if (this->data().hasProperty(map["connectivity"]))
    {
        auto c = this->data().property(map["connectivity"]).asA<Connectivity>();

        auto int_atoms0 = _to_int(atoms0);
        auto int_atoms1 = _to_int(atoms1);

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto angles01 = _filter(angles, int_atoms0, 0);
        angles01 = _filter(angles01, int_atoms1, 1);

        auto angles12 = _filter(angles, int_atoms0, 1);
        angles12 = _filter(angles12, int_atoms1, 2);

        auto angles21 = _filter(angles, int_atoms0, 2);
        angles21 = _filter(angles21, int_atoms1, 1);

        auto angles10 = _filter(angles, int_atoms0, 1);
        angles10 = _filter(angles10, int_atoms1, 0);

        angs = _from_int_set(angles01 + angles12 + angles21 + angles10);
    }
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
            "You can only create an angle from atoms in the same molecule. "
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

        auto angles = _to_int_set(c.getAngles(), this->data().info());

        auto angles012 = _filter(angles, int_atoms0, 0);
        angles012 = _filter(angles012, int_atoms1, 1);
        angles012 = _filter(angles012, int_atoms2, 2);

        auto angles210 = _filter(angles, int_atoms0, 2);
        angles210 = _filter(angles210, int_atoms1, 1);
        angles210 = _filter(angles210, int_atoms2, 0);

        angs = _from_int_set(angles012 + angles210);
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

QList<SireUnits::Dimension::Angle> SelectorAngle::measures(const PropertyMap &map) const
{
    return this->sizes(map);
}

QList<SireUnits::Dimension::Angle> SelectorAngle::measures() const
{
    return this->sizes();
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
