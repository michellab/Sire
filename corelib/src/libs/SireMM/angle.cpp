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

#include "angle.h"
#include "selectorangle.h"

#include "threeatomfunctions.h"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"

#include "SireCAS/symbol.h"
#include "SireCAS/values.h"
#include "SireCAS/expression.h"

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireStream;
using namespace SireError;
using namespace SireUnits;

static const RegisterMetaType<Angle> r_angle;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const Angle &angle)
{
    writeHeader(ds, r_angle, 1);
    ds << angle.ang << static_cast<const MoleculeView&>(angle);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, Angle &angle)
{
    auto v = readHeader(ds, r_angle);

    if (v == 1)
    {
        ds >> angle.ang >> static_cast<MoleculeView&>(angle);
    }
    else
        throw SireStream::version_error(v, "1", r_angle, CODELOC);

    return ds;
}

Angle::Angle() : ConcreteProperty<Angle, MoleculeView>()
{}


Angle::Angle(const Atom &atom0, const Atom &atom1, const Atom &atom2)
     : ConcreteProperty<Angle, MoleculeView>()
{
    if ((not atom0.isSameMolecule(atom1)) or (not atom0.isSameMolecule(atom2)))
    {
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Angle from three atoms in the same molecule. "
            "%1-%2-%3 are from different molecules (%4-%5-%6)")
                .arg(atom0.toString()).arg(atom1.toString())
                .arg(atom2.toString())
                .arg(atom0.molecule().toString())
                .arg(atom1.molecule().toString())
                .arg(atom2.molecule().toString()), CODELOC);
    }

    this->operator=(Angle(atom0.data(),
                          AngleID(atom0.index(),
                                  atom1.index(),
                                  atom2.index())));
}

Angle::Angle(const MoleculeView &molview,
             const AtomID &atom0, const AtomID &atom1,
             const AtomID &atom2)
     : ConcreteProperty<Angle, MoleculeView>()
{
    this->operator=(Angle(molview.atom(atom0),
                          molview.atom(atom1),
                          molview.atom(atom2)));
}

Angle::Angle(const MoleculeData &moldata,
             const AtomID &atm0, const AtomID &atm1,
             const AtomID &atm2)
     : ConcreteProperty<Angle, MoleculeView>()
{
    this->operator=(Angle(Atom(moldata, atm0),
                          Atom(moldata, atm1),
                          Atom(moldata, atm2)));
}

Angle::Angle(const MoleculeData &moldata, const AngleID &angle)
      : ConcreteProperty<Angle, MoleculeView>(moldata)
{
    auto atomidx0 = moldata.info().atomIdx(angle.atom0());
    auto atomidx1 = moldata.info().atomIdx(angle.atom1());
    auto atomidx2 = moldata.info().atomIdx(angle.atom2());

    if (atomidx0 > atomidx2)
    {
        qSwap(atomidx0, atomidx2);
    }
    else if ((atomidx0 == atomidx1) or
             (atomidx0 == atomidx2) or
             (atomidx1 == atomidx2))
    {
        throw SireMol::duplicate_atom(QObject::tr(
            "You cannot make a Angle out of identical atoms. %1-%2-%3")
                .arg(atomidx0.toString())
                .arg(atomidx1.toString())
                .arg(atomidx2.toString()), CODELOC);
    }

    ang = AngleID(atomidx0, atomidx1, atomidx2);
}

Angle::Angle(const Angle &other)
     : ConcreteProperty<Angle, MoleculeView>(other), ang(other.ang)
{}

Angle::~Angle()
{}

const char* Angle::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Angle>());
}

Angle& Angle::operator=(const Angle &other)
{
    if (this != &other)
    {
        ang = other.ang;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool Angle::operator==(const Angle &other) const
{
    return ang == other.ang and MoleculeView::operator==(other);
}

bool Angle::operator!=(const Angle &other) const
{
    return not operator==(other);
}

MolViewPtr Angle::toSelector() const
{
    return SelectorAngle(*this);
}

QString Angle::toString() const
{
    if (this->isNull())
        return QObject::tr("Angle::null");

    auto a0 = this->atom0();
    auto a1 = this->atom1();
    auto a2 = this->atom2();

    return QObject::tr("Angle( %1:%2 <= %3:%4 => %5:%6 )")
            .arg(a0.name()).arg(a0.number())
            .arg(a1.name()).arg(a1.number())
            .arg(a2.name()).arg(a2.number());
}

Atom Angle::atom0() const
{
    return Atom(*this, ang.atom0());
}

Atom Angle::atom1() const
{
    return Atom(*this, ang.atom1());
}

Atom Angle::atom2() const
{
    return Atom(*this, ang.atom2());
}

Evaluator Angle::evaluate() const
{
    return Evaluator(*this);
}

Mover<Angle> Angle::move() const
{
    return Mover<Angle>(*this);
}

AngleID Angle::ID() const
{
    return ang;
}

bool Angle::isEmpty() const
{
    return this->isNull();
}

bool Angle::selectedAll() const
{
    return this->data().info().nAtoms() == 3;
}

AtomSelection Angle::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    s = s.select(ang.atom0());
    s = s.select(ang.atom1());
    s = s.select(ang.atom2());

    return s;
}

bool Angle::hasConnectivity() const
{
    return this->data().hasProperty("connectivity");
}

const Connectivity& Angle::getConnectivity() const
{
    return this->data().property("connectivity").asA<Connectivity>();
}

bool Angle::hasProperty(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().hasProperty(ang, key);
    }

    return false;
}

bool Angle::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool Angle::hasMetadata(const PropertyName &key,
                        const PropertyName &metakey) const
{
    return false;
}

QStringList Angle::propertyKeys() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().propertyKeys(ang);
    }

    return QStringList();
}

QStringList Angle::metadataKeys() const
{
    return QStringList();
}

QStringList Angle::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

Properties Angle::properties() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().properties(ang);
    }

    return Properties();
}

const Property& Angle::property(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(ang, key);
    }

    if (key.hasValue())
        return key.value();

    throw SireBase::missing_property(QObject::tr(
        "Angle %1 has no property at key %2.")
            .arg(this->toString()).arg(key.source()),
                CODELOC);

    return key.value();
}

const Property& Angle::property(const PropertyName &key,
                                const Property &default_value) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(ang, key, default_value);
    }

    if (key.hasValue())
        return key.value();

    return default_value;
}

SireUnits::Dimension::Angle Angle::size(const PropertyMap &map) const
{
    return ang.size(*this, map);
}

SireUnits::Dimension::Angle Angle::size() const
{
    return this->size(PropertyMap());
}

SireCAS::Expression Angle::potential() const
{
    return this->potential(PropertyMap());
}

Expression Angle::potential(const PropertyMap &map) const
{
    try
    {
        auto funcs = this->data().property(map["angle"]).asA<ThreeAtomFunctions>();
        return funcs.potential(ang);
    }
    catch(const SireError::exception&)
    {}

    return Expression();
}

SireUnits::Dimension::MolarEnergy Angle::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy Angle::energy(const PropertyMap &map) const
{
    auto pot = this->potential(map);
    auto s = this->size(map);

    Values vals(Symbol("theta")==s.to(radians));
    return pot.evaluate(vals) * kcal_per_mol;
}

namespace SireMol
{
    template class Mover<Angle>;
}
