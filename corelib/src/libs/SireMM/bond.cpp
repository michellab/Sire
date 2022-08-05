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

#include "bond.h"
#include "selectorbond.h"

#include "twoatomfunctions.h"

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
using namespace SireUnits::Dimension;

static const RegisterMetaType<Bond> r_bond;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const Bond &bond)
{
    writeHeader(ds, r_bond, 1);
    ds << bond.bnd << static_cast<const MoleculeView&>(bond);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, Bond &bond)
{
    auto v = readHeader(ds, r_bond);

    if (v == 1)
    {
        ds >> bond.bnd >> static_cast<MoleculeView&>(bond);
    }
    else
        throw SireStream::version_error(v, "1", r_bond, CODELOC);

    return ds;
}

Bond::Bond() : ConcreteProperty<Bond, MoleculeView>()
{}


Bond::Bond(const Atom &atom0, const Atom &atom1)
     : ConcreteProperty<Bond, MoleculeView>()
{
    if (not atom0.isSameMolecule(atom1))
    {
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a bond from two atoms in the same molecule. "
            "%1 and %2 are from different molecules (%3 and %4)")
                .arg(atom0.toString()).arg(atom1.toString())
                .arg(atom0.molecule().toString())
                .arg(atom1.molecule().toString()), CODELOC);
    }

    this->operator=(Bond(atom0.data(),
                         BondID(atom0.index(), atom1.index())));
}

Bond::Bond(const MoleculeView &molview,
           const AtomID &atom0, const AtomID &atom1)
     : ConcreteProperty<Bond, MoleculeView>(molview)
{
    this->operator=(Bond(molview.atom(atom0), molview.atom(atom1)));
}

Bond::Bond(const MoleculeData &moldata,
           const AtomID &atm0, const AtomID &atm1)
     : ConcreteProperty<Bond, MoleculeView>()
{
    this->operator=(Bond(Atom(moldata, atm0), Atom(moldata, atm1)));
}

Bond::Bond(const MoleculeData &moldata, const BondID &bond)
     : ConcreteProperty<Bond, MoleculeView>(moldata)
{
    auto atomidx0 = moldata.info().atomIdx(bond.atom0());
    auto atomidx1 = moldata.info().atomIdx(bond.atom1());

    if (atomidx0 > atomidx1)
    {
        qSwap(atomidx0, atomidx1);
    }
    else if (atomidx0 == atomidx1)
    {
        throw SireMol::duplicate_atom(QObject::tr(
            "You cannot make a bond out of two identical atoms. %1")
                .arg(atomidx0.toString()), CODELOC);
    }

    bnd = BondID(atomidx0, atomidx1);
}

Bond::Bond(const Bond &other)
     : ConcreteProperty<Bond, MoleculeView>(other), bnd(other.bnd)
{}

Bond::~Bond()
{}

const char* Bond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Bond>());
}

Bond& Bond::operator=(const Bond &other)
{
    if (this != &other)
    {
        bnd = other.bnd;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool Bond::operator==(const Bond &other) const
{
    return bnd == other.bnd and MoleculeView::operator==(other);
}

bool Bond::operator!=(const Bond &other) const
{
    return not operator==(other);
}

MolViewPtr Bond::toSelector() const
{
    return SelectorBond(*this);
}

QString Bond::toString() const
{
    if (this->isNull())
        return QObject::tr("Bond::null");

    auto a0 = this->atom0();
    auto a1 = this->atom1();

    return QObject::tr("Bond( %1:%2 => %3:%4 )")
            .arg(a0.name()).arg(a0.number())
            .arg(a1.name()).arg(a1.number());
}

Atom Bond::atom0() const
{
    return Atom(*this, bnd.atom0());
}

Atom Bond::atom1() const
{
    return Atom(*this, bnd.atom1());
}

Evaluator Bond::evaluate() const
{
    return Evaluator(*this);
}

Mover<Bond> Bond::move() const
{
    return Mover<Bond>(*this);
}

BondID Bond::ID() const
{
    return bnd;
}

bool Bond::isEmpty() const
{
    return this->isNull();
}

bool Bond::selectedAll() const
{
    return this->data().info().nAtoms() == 2;
}

AtomSelection Bond::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    s = s.select(bnd.atom0());
    s = s.select(bnd.atom1());

    return s;
}

bool Bond::hasConnectivity() const
{
    return this->data().hasProperty("connectivity");
}

const Connectivity& Bond::getConnectivity() const
{
    return this->data().property("connectivity").asA<Connectivity>();
}

bool Bond::hasProperty(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().hasProperty(bnd, key);
    }

    return false;
}

bool Bond::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool Bond::hasMetadata(const PropertyName &key,
                       const PropertyName &metakey) const
{
    return false;
}

QStringList Bond::propertyKeys() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().propertyKeys(bnd);
    }

    return QStringList();
}

QStringList Bond::metadataKeys() const
{
    return QStringList();
}

QStringList Bond::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

Properties Bond::properties() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().properties(bnd);
    }

    return Properties();
}

const Property& Bond::property(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(bnd, key);
    }

    if (key.hasValue())
        return key.value();

    throw SireBase::missing_property(QObject::tr(
        "Bond %1 has no property at key %2.")
            .arg(this->toString()).arg(key.source()),
                CODELOC);

    return key.value();
}

const Property& Bond::property(const PropertyName &key,
                               const Property &default_value) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(bnd, key, default_value);
    }

    if (key.hasValue())
        return key.value();

    return default_value;
}

SireUnits::Dimension::Length Bond::length(const PropertyMap &map) const
{
    return bnd.length(*this, map) * angstrom;
}

SireUnits::Dimension::Length Bond::length() const
{
    return this->length(PropertyMap());
}

SireCAS::Expression Bond::potential() const
{
    return this->potential(PropertyMap());
}

Expression Bond::potential(const PropertyMap &map) const
{
    try
    {
        auto funcs = this->data().property(map["bond"]).asA<TwoAtomFunctions>();
        return funcs.potential(bnd);
    }
    catch(const SireError::exception&)
    {}

    return Expression();
}

SireUnits::Dimension::MolarEnergy Bond::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy Bond::energy(const PropertyMap &map) const
{
    auto pot = this->potential(map);
    auto l = this->length(map);

    Values vals(Symbol("r")==l.to(angstrom));
    return pot.evaluate(vals) * kcal_per_mol;
}

namespace SireMol
{
    template class Mover<Bond>;
}
