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

#include "improper.h"
#include "selectorimproper.h"

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

static const RegisterMetaType<Improper> r_improper;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const Improper &improper)
{
    writeHeader(ds, r_improper, 1);
    ds << improper.imp << static_cast<const MoleculeView&>(improper);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, Improper &improper)
{
    auto v = readHeader(ds, r_improper);

    if (v == 1)
    {
        ds >> improper.imp >> static_cast<MoleculeView&>(improper);
    }
    else
        throw SireStream::version_error(v, "1", r_improper, CODELOC);

    return ds;
}

Improper::Improper() : ConcreteProperty<Improper, MoleculeView>()
{}


Improper::Improper(const Atom &atom0, const Atom &atom1,
                   const Atom &atom2, const Atom &atom3)
         : ConcreteProperty<Improper, MoleculeView>()
{
    if ((not atom0.isSameMolecule(atom1)) or
        (not atom0.isSameMolecule(atom2)) or
        (not atom0.isSameMolecule(atom3)))
    {
        throw SireError::incompatible_error(QObject::tr(
            "You can only create an Improper from four atoms in the same molecule. "
            "%1-%2-%3-%4 are from different molecules (%5-%6-%7-%8)")
                .arg(atom0.toString()).arg(atom1.toString())
                .arg(atom2.toString()).arg(atom3.toString())
                .arg(atom0.molecule().toString())
                .arg(atom1.molecule().toString())
                .arg(atom2.molecule().toString())
                .arg(atom3.molecule().toString()), CODELOC);
    }

    this->operator=(Improper(atom0.data(),
                          ImproperID(atom0.index(),
                                     atom1.index(),
                                     atom2.index(),
                                     atom3.index())));
}

Improper::Improper(const MoleculeView &molview,
                   const AtomID &atom0, const AtomID &atom1,
                   const AtomID &atom2, const AtomID &atom3)
     : ConcreteProperty<Improper, MoleculeView>()
{
    this->operator=(Improper(molview.atom(atom0),
                             molview.atom(atom1),
                             molview.atom(atom2),
                             molview.atom(atom3)));
}

Improper::Improper(const MoleculeData &moldata,
                   const AtomID &atm0, const AtomID &atm1,
                   const AtomID &atm2, const AtomID &atm3)
     : ConcreteProperty<Improper, MoleculeView>()
{
    this->operator=(Improper(Atom(moldata, atm0),
                             Atom(moldata, atm1),
                             Atom(moldata, atm2),
                             Atom(moldata, atm3)));
}

Improper::Improper(const MoleculeData &moldata, const ImproperID &improper)
         : ConcreteProperty<Improper, MoleculeView>(moldata)
{
    auto atomidx0 = moldata.info().atomIdx(improper.atom0());
    auto atomidx1 = moldata.info().atomIdx(improper.atom1());
    auto atomidx2 = moldata.info().atomIdx(improper.atom2());
    auto atomidx3 = moldata.info().atomIdx(improper.atom3());

    if ((atomidx0 == atomidx1) or
        (atomidx0 == atomidx2) or
        (atomidx0 == atomidx3) or
        (atomidx1 == atomidx2) or
        (atomidx1 == atomidx3) or
        (atomidx2 == atomidx3))
    {
        throw SireMol::duplicate_atom(QObject::tr(
            "You cannot make an Improper out of identical atoms. %1-%2-%3-%4")
                .arg(atomidx0.toString())
                .arg(atomidx1.toString())
                .arg(atomidx2.toString())
                .arg(atomidx3.toString()), CODELOC);
    }

    imp = ImproperID(atomidx0, atomidx1, atomidx2, atomidx3);
}

Improper::Improper(const Improper &other)
         : ConcreteProperty<Improper, MoleculeView>(other), imp(other.imp)
{}

Improper::~Improper()
{}

const char* Improper::typeName()
{
    return QMetaType::typeName(qMetaTypeId<Improper>());
}

Improper& Improper::operator=(const Improper &other)
{
    if (this != &other)
    {
        imp = other.imp;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool Improper::operator==(const Improper &other) const
{
    return imp == other.imp and MoleculeView::operator==(other);
}

bool Improper::operator!=(const Improper &other) const
{
    return not operator==(other);
}

MolViewPtr Improper::toSelector() const
{
    return SelectorImproper(*this);
}

QString Improper::toString() const
{
    if (this->isNull())
        return QObject::tr("Improper::null");

    auto a0 = this->atom0();
    auto a1 = this->atom1();
    auto a2 = this->atom2();
    auto a3 = this->atom3();

    return QObject::tr("Improper( %1:%2 >= %3:%4 = %5:%6 <= %7:%8 )")
            .arg(a0.name()).arg(a0.number())
            .arg(a1.name()).arg(a1.number())
            .arg(a2.name()).arg(a2.number())
            .arg(a3.name()).arg(a3.number());
}

Atom Improper::atom0() const
{
    return Atom(*this, imp.atom0());
}

Atom Improper::atom1() const
{
    return Atom(*this, imp.atom1());
}

Atom Improper::atom2() const
{
    return Atom(*this, imp.atom2());
}

Atom Improper::atom3() const
{
    return Atom(*this, imp.atom3());
}

Evaluator Improper::evaluate() const
{
    return Evaluator(*this);
}

Mover<Improper> Improper::move() const
{
    return Mover<Improper>(*this);
}

ImproperID Improper::ID() const
{
    return imp;
}

bool Improper::isEmpty() const
{
    return this->isNull();
}

bool Improper::selectedAll() const
{
    return this->data().info().nAtoms() == 4;
}

AtomSelection Improper::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    s = s.select(imp.atom0());
    s = s.select(imp.atom1());
    s = s.select(imp.atom2());
    s = s.select(imp.atom3());

    return s;
}

bool Improper::hasConnectivity() const
{
    return this->data().hasProperty("connectivity");
}

const Connectivity& Improper::getConnectivity() const
{
    return this->data().property("connectivity").asA<Connectivity>();
}

bool Improper::hasProperty(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().hasProperty(imp, key);
    }

    return false;
}

bool Improper::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool Improper::hasMetadata(const PropertyName &key,
                           const PropertyName &metakey) const
{
    return false;
}

QStringList Improper::propertyKeys() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().propertyKeys(imp);
    }

    return QStringList();
}

QStringList Improper::metadataKeys() const
{
    return QStringList();
}

QStringList Improper::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

Properties Improper::properties() const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().properties(imp);
    }

    return Properties();
}

const Property& Improper::property(const PropertyName &key) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(imp, key);
    }

    if (key.hasValue())
        return key.value();

    throw SireBase::missing_property(QObject::tr(
        "Improper %1 has no property at key %2.")
            .arg(this->toString()).arg(key.source()),
                CODELOC);

    return key.value();
}

const Property& Improper::property(const PropertyName &key,
                                   const Property &default_value) const
{
    if (this->hasConnectivity())
    {
        return this->getConnectivity().property(imp, key, default_value);
    }

    if (key.hasValue())
        return key.value();

    return default_value;
}

SireUnits::Dimension::Angle Improper::size(const PropertyMap &map) const
{
    return imp.size(*this, map);
}

SireUnits::Dimension::Angle Improper::size() const
{
    return this->size(PropertyMap());
}

SireCAS::Expression Improper::potential() const
{
    return this->potential(PropertyMap());
}

Expression Improper::potential(const PropertyMap &map) const
{
    try
    {
        auto funcs = this->data().property(map["improper"]).asA<FourAtomFunctions>();
        return funcs.potential(imp);
    }
    catch(const SireError::exception&)
    {}

    return Expression();
}

SireUnits::Dimension::MolarEnergy Improper::energy() const
{
    return this->energy(PropertyMap());
}

SireUnits::Dimension::MolarEnergy Improper::energy(const PropertyMap &map) const
{
    auto pot = this->potential(map);
    auto s = this->size(map);

    Values vals(Symbol("phi")==s.to(radians), Symbol("theta")==s.to(radians));
    return pot.evaluate(vals) * kcal_per_mol;
}

namespace SireMol
{
    template class Mover<Improper>;
}
