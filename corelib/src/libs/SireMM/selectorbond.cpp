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

#include "selectorbond.h"

#include "twoatomfunctions.h"

#include "SireMol/molecule.h"
#include "SireMol/mover.hpp"
#include "SireMol/selector.hpp"

#include "SireCAS/symbol.h"
#include "SireCAS/values.h"

#include "SireID/index.h"

#include "SireBase/errors.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireCAS;
using namespace SireID;
using namespace SireStream;
using namespace SireError;
using namespace SireUnits;
using namespace SireUnits::Dimension;

static const RegisterMetaType<SelectorBond> r_sbond;

SIREMM_EXPORT QDataStream& operator<<(QDataStream &ds, const SelectorBond &sbond)
{
    writeHeader(ds, r_sbond, 1);
    SharedDataStream sds(ds);
    sds << sbond.bnds << static_cast<const MoleculeView&>(sbond);
    return ds;
}

SIREMM_EXPORT QDataStream& operator>>(QDataStream &ds, SelectorBond &sbond)
{
    auto v = readHeader(ds, r_sbond);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> sbond.bnds >> static_cast<MoleculeView&>(sbond);
    }
    else
        throw SireStream::version_error(v, "1", r_sbond, CODELOC);

    return ds;
}

SelectorBond::SelectorBond() : ConcreteProperty<SelectorBond, MoleculeView>()
{}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    if (mol.data().hasProperty(map["connectivity"]))
    {
        auto c = mol.data().property(map["connectivity"]).asA<Connectivity>();

        auto bonds = c.getBonds();

        if (mol.selectedAll())
        {
            bnds = bonds;
        }
        else
        {
            const auto s = mol.selection();

            for (const auto &bond : bonds)
            {
                if (s.selected(bond.atom0()) and s.selected(bond.atom1()))
                {
                    bnds.append(bond);
                }
            }
        }
    }
}

SelectorBond::SelectorBond(const MoleculeView &mol,
                           const QList<BondID> &bonds)
             : ConcreteProperty<SelectorBond, MoleculeView>(mol)
{
    const auto s = mol.selection();

    for (const auto &bond : bonds)
    {
        BondID b(mol.data().info().atomIdx(bond.atom0()),
                 mol.data().info().atomIdx(bond.atom1()));

        if (s.selected(b.atom0()) and s.selected(b.atom1()))
        {
            bnds.append(b);
        }
    }
}

SelectorBond::SelectorBond(const MoleculeData &moldata,
                           const SireBase::PropertyMap &map)
     : ConcreteProperty<SelectorBond, MoleculeView>()
{
    this->operator=(SelectorBond(Molecule(moldata), map));
}

SelectorBond::SelectorBond(const SelectorBond &other)
     : ConcreteProperty<SelectorBond, MoleculeView>(other), bnds(other.bnds)
{}

SelectorBond::~SelectorBond()
{}

const char* SelectorBond::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SelectorBond>());
}

SelectorBond& SelectorBond::operator=(const SelectorBond &other)
{
    if (this != &other)
    {
        bnds = other.bnds;
        MoleculeView::operator=(other);
    }

    return *this;
}

bool SelectorBond::operator==(const SelectorBond &other) const
{
    return bnds == other.bnds and MoleculeView::operator==(other);
}

bool SelectorBond::operator!=(const SelectorBond &other) const
{
    return not operator==(other);
}

int SelectorBond::count() const
{
    return this->bnds.count();
}

int SelectorBond::size() const
{
    return this->count();
}

int SelectorBond::nViews() const
{
    return this->count();
}

QString SelectorBond::toString() const
{
    if (this->isNull())
        return QObject::tr("SelectorBond::null");

    QStringList parts;

    if (this->count() <= 10)
    {
        for (int i=0; i<10; ++i)
        {
            parts.append(this->operator()(i).toString());
        }
    }
    else
    {
        for (int i=0; i<5; ++i)
        {
            parts.append(this->operator()(i).toString());
        }

        parts.append("...");

        for (int i=this->count()-5; i<this->count(); ++i)
        {
            parts.append(this->operator()(i).toString());
        }
    }

    return QObject::tr("SelectorMol( size=%1\n%2\n)")
                .arg(this->count()).arg(parts.join("\n"));
}

MolViewPtr SelectorBond::operator[](int i) const
{
    return this->operator()(i);
}

Bond SelectorBond::operator()(int i) const
{
    auto bnd = bnds.at(Index(i).map(bnds.count()));
    return Bond(this->data(), bnd);
}

SelectorBond SelectorBond::operator()(int i, int j) const
{
    i = Index(i).map(bnds.count());
    j = Index(j).map(bnds.count());

    SelectorBond ret(*this);
    ret.bnds.clear();

    if (i <= j)
    {
        for ( ; i<=j; ++i)
        {
            ret.bnds.append(this->bnds.at(i));
        }
    }
    else
    {
        for ( ; i >= j; --i)
        {
            ret.bnds.append(this->bnds.at(i));
        }
    }

    return ret;
}

QList<BondID> SelectorBond::IDs() const
{
    return bnds;
}

bool SelectorBond::isEmpty() const
{
    return this->bnds.isEmpty();
}

bool SelectorBond::selectedAll() const
{
    return this->selection().selectedAll();
}

AtomSelection SelectorBond::selection() const
{
    if (this->isNull())
        return AtomSelection();

    auto s = AtomSelection(this->data());
    s = s.deselectAll();

    for (const auto &bnd : bnds)
    {
        s = s.select(bnd.atom0());
        s = s.select(bnd.atom1());
    }

    return s;
}

bool SelectorBond::hasProperty(const PropertyName &key) const
{
    for (int i=0; i<this->count(); ++i)
    {
        if (this->operator()(i).hasProperty(key))
            return true;
    }

    return false;
}

bool SelectorBond::hasMetadata(const PropertyName &key) const
{
    return false;
}

bool SelectorBond::hasMetadata(const PropertyName &key,
                               const PropertyName &metakey) const
{
    return false;
}

QStringList SelectorBond::propertyKeys() const
{
    QSet<QString> keys;

    for (int i=0; i<this->count(); ++i)
    {
        for (const auto &k : this->operator()(i).propertyKeys())
        {
            keys.insert(k);
        }
    }

    QStringList ret(keys.constBegin(), keys.constEnd());
    return ret;
}

QStringList SelectorBond::metadataKeys() const
{
    return QStringList();
}

QStringList SelectorBond::metadataKeys(const PropertyName &key) const
{
    return QStringList();
}

QList<Properties> SelectorBond::properties() const
{
    QList<Properties> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).properties());
    }

    return props;
}

QList<PropertyPtr> SelectorBond::property(const PropertyName &key) const
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
            "None of the bonds in this container have a property called %1.")
                .arg(key.source()), CODELOC);

    return props;
}

QList<PropertyPtr> SelectorBond::property(const PropertyName &key,
                                          const Property &default_value) const
{
    QList<PropertyPtr> props;

    for (int i=0; i<this->count(); ++i)
    {
        props.append(this->operator()(i).property(key, default_value));
    }

    return props;
}

QList<Length> SelectorBond::lengths(const PropertyMap &map) const
{
    QList<Length> l;

    for (int i=0; i<this->count(); ++i)
    {
        l.append(this->operator()(i).length(map));
    }

    return l;
}

QList<Length> SelectorBond::lengths() const
{
    return this->lengths(PropertyMap());
}

QList<Expression> SelectorBond::potentials() const
{
    return this->potentials(PropertyMap());
}

QList<Expression> SelectorBond::potentials(const PropertyMap &map) const
{
    QList<Expression> p;

    for (int i=0; i<this->count(); ++i)
    {
        p.append(this->operator()(i).potential(map));
    }

    return p;
}

QList<MolarEnergy> SelectorBond::energies() const
{
    return this->energies(PropertyMap());
}

QList<MolarEnergy> SelectorBond::energies(const PropertyMap &map) const
{
    QList<MolarEnergy> nrgs;

    for (int i=0; i<this->count(); ++i)
    {
        nrgs.append(this->operator()(i).energy(map));
    }

    return nrgs;
}

MolarEnergy SelectorBond::energy() const
{
    return this->energy(PropertyMap());
}

MolarEnergy SelectorBond::energy(const PropertyMap &map) const
{
    MolarEnergy nrg(0);

    for (int i=0; i<this->count(); ++i)
    {
        nrg += this->operator()(i).energy(map);
    }

    return nrg;
}
