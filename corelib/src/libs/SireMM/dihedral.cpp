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

#include "dihedral.h"
#include "selectordihedral.h"

#include "fouratomfunctions.h"

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

static const RegisterMetaType<Dihedral> r_dihedral;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const Dihedral &dihedral)
{
    writeHeader(ds, r_dihedral, 1);
    ds << dihedral.dih << static_cast<const MoleculeView&>(dihedral);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, Dihedral &dihedral)
{
    auto v = readHeader(ds, r_dihedral);

    if (v == 1)
    {
        ds >> dihedral.dih >> static_cast<MoleculeView&>(dihedral);
    }
    else
        throw SireStream::version_error(v, "1", r_dihedral, CODELOC);

    return ds;
}

Dihedral::Dihedral() : ConcreteProperty<Dihedral, MoleculeView>()
{}


Dihedral::Dihedral(const Atom &atom0, const Atom &atom1,
                   const Atom &atom2, const Atom &atom3)
     : ConcreteProperty<Dihedral, MoleculeView>()
{
    if ((not atom0.isSameMolecule(atom1)) or
        (not atom0.isSameMolecule(atom2)) or
        (not atom0.isSameMolecule(atom3)))
    {
        throw SireError::incompatible_error(QObject::tr(
            "You can only create a Dihedral from four atoms in the same molecule. "
            "%1-%2-%3-%4 are from different molecules (%5-%6-%7-%8)")
                .arg(atom0.toString()).arg(atom1.toString())
                .arg(atom2.toString()).arg(atom3.toString())
                .arg(atom0.molecule().toString())
                .arg(atom1.molecule().toString())
                .arg(atom2.molecule().toString())
                .arg(atom3.molecule().toString()), CODELOC);
    }

    this->operator=(Dihedral(atom0.data(),
                          DihedralID(atom0.index(),
                                     atom1.index(),
                                     atom2.index(),
                                     atom3.index())));
}

Dihedral::Dihedral(const MoleculeView &molview,
                   const AtomID &atom0, const AtomID &atom1,
                   const AtomID &atom2, const AtomID &atom3)
     : ConcreteProperty<Dihedral, MoleculeView>()
{
    this->operator=(Dihedral(molview.atom(atom0),
                             molview.atom(atom1),
                             molview.atom(atom2),
                             molview.atom(atom3)));
}

Dihedral::Dihedral(const MoleculeData &moldata,
                   const AtomID &atm0, const AtomID &atm1,
                   const AtomID &atm2, const AtomID &atm3)
     : ConcreteProperty<Dihedral, MoleculeView>()
{
    this->operator=(Dihedral(Atom(moldata, atm0),
                             Atom(moldata, atm1),
                             Atom(moldata, atm2),
                             Atom(moldata, atm3)));
}

Dihedral::Dihedral(const MoleculeData &moldata, const DihedralID &dihedral)
         : ConcreteProperty<Dihedral, MoleculeView>(moldata)
{
    auto atomidx0 = moldata.info().atomIdx(dihedral.atom0());
    auto atomidx1 = moldata.info().atomIdx(dihedral.atom1());
    auto atomidx2 = moldata.info().atomIdx(dihedral.atom2());
    auto atomidx3 = moldata.info().atomIdx(dihedral.atom3());

    if (atomidx0 > atomidx3)
    {
        qSwap(atomidx0, atomidx3);
        qSwap(atomidx1, atomidx2);
    }
    else if ((atomidx0 == atomidx1) or
             (atomidx1 == atomidx2) or
             (atomidx2 == atomidx3))
    {
        // note that atomidx0 == atomidx3 in a ring
        throw SireMol::duplicate_atom(QObject::tr(
            "You cannot make a Dihedral out of identical atoms. %1-%2-%3-%4")
                .arg(atomidx0.toString())
                .arg(atomidx1.toString())
                .arg(atomidx2.toString())
                .arg(atomidx3.toString()), CODELOC);
    }

    dih = DihedralID(atomidx0, atomidx1, atomidx2, atomidx3);
}

Dihedral::Dihedral(const Dihedral &other)
         : ConcreteProperty<Dihedral, MoleculeView>(other), dih(other.dih)
{}

Dihedral::~Dihedral()
{}

const char* Dihedral::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Dihedral>());
}

Dihedral& Dihedral::operator=(const Dihedral &other)
{
    if (this != &other)
    {
        dih = other.dih;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool Dihedral::operator==(const Dihedral &other) const
{
    return dih == other.dih and MoleculeView::operator==(other);
}

bool Dihedral::operator!=(const Dihedral &other) const
{
    return not operator==(other);
}

MolViewPtr Dihedral::toSelector() const
{
    return SelectorDihedral(*this);
}

QString Dihedral::toString() const
{
    if (this->isNull())
        return QObject::tr("Dihedral::null");

    auto a0 = this->atom0();
    auto a1 = this->atom1();
    auto a2 = this->atom2();
    auto a3 = this->atom3();

    return QObject::tr("Dihedral( %1:%2 <= %3:%4 = %5:%6 => %7:%8 )")
            .arg(a0.name()).arg(a0.number())
            .arg(a1.name()).arg(a1.number())
            .arg(a2.name()).arg(a2.number())
            .arg(a3.name()).arg(a3.number());
}

Atom Dihedral::atom0() const
{
    return Atom(*this, dih.atom0());
}

Atom Dihedral::atom1() const
{
    return Atom(*this, dih.atom1());
}

Atom Dihedral::atom2() const
{
    return Atom(*this, dih.atom2());
}

Atom Dihedral::atom3() const
{
    return Atom(*this, dih.atom3());
}

Evaluator Dihedral::evaluate() const
{
    return Evaluator(*this);
}

Mover<Dihedral> Dihedral::move() const
{
    return Mover<Dihedral>(*this);
}

DihedralID Dihedral::ID() const
{
    return dih;
}

bool Dihedral::isEmpty() const
{
    return this->isNull();
}

bool Dihedral::selectedAll() const
{
    return this->data().info().nAtoms() == 4;
}

AtomSelection Dihedral::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    s = s.select(dih.atom0());
    s = s.select(dih.atom1());
    s = s.select(dih.atom2());
    s = s.select(dih.atom3());

    return s;
}

bool Dihedral::hasConnectivity() const
{
    return this->data().hasProperty("connectivity");
}

const Connectivity& Dihedral::getConnectivity() const
{
    return this->data().property("connectivity").asA<Connectivity>();
}

bool Dihedral::hasProperty(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().hasProperty(dih, key);
    }

    return false;
}

bool Dihedral::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool Dihedral::hasMetadata(const PropertyName &key,
                           const PropertyName &metakey) const
{
    return false;
}

QStringList Dihedral::propertyKeys() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().propertyKeys(dih);
    }

    return QStringList();
}

QStringList Dihedral::metadataKeys() const
{
    return QStringList();
}

QStringList Dihedral::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

Properties Dihedral::properties() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().properties(dih);
    }

    return Properties();
}

const Property& Dihedral::property(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(dih, key);
    }

    if (key.hasValue())
        return key.value();

    throw SireBase::missing_property(QObject::tr(
        "Dihedral %1 has no property at key %2.")
            .arg(this->toString()).arg(key.source()),
                CODELOC);

    return key.value();
}

const Property& Dihedral::property(const PropertyName &key,
                                   const Property &default_value) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(dih, key, default_value);
    }

    if (key.hasValue())
        return key.value();

    return default_value;
}

SireUnits::Dimension::Angle Dihedral::size(const PropertyMap &map) const
{
    return dih.size(*this, map);
}

SireUnits::Dimension::Angle Dihedral::size() const
{
    return this->size(PropertyMap());
}

SireUnits::Dimension::Angle Dihedral::measure(const PropertyMap &map) const
{
    return this->size(map);
}

SireUnits::Dimension::Angle Dihedral::measure() const
{
    return this->size();
}

SireCAS::Expression Dihedral::potential() const
{
    return this->potential(PropertyMap());
}

Expression Dihedral::potential(const PropertyMap &map) const
{
    try
    {
        auto funcs = this->data().property(map["dihedral"]).asA<FourAtomFunctions>();
        return funcs.potential(dih);
    }
    catch(const SireError::exception&)
    {}

    return Expression();
}

SireUnits::Dimension::MolarEnergy Dihedral::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy Dihedral::energy(const PropertyMap &map) const
{
    auto pot = this->potential(map);
    auto s = this->size(map);

    Values vals(Symbol("phi")==s.to(radians));
    return pot.evaluate(vals) * kcal_per_mol;
}

namespace SireMol
{
    template class Mover<Dihedral>;
}
