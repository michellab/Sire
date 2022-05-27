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

#ifndef SIREMM_BOND_H
#define SIREMM_BOND_H

#include "SireMol/atom.h"
#include "SireMol/bondid.h"
#include "SireMol/angleid.h"
#include "SireMol/dihedralid.h"
#include "SireMol/improperid.h"
#include "SireMol/connectivity.h"
#include "SireMol/mover.hpp"
#include "SireMol/evaluator.h"

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class Bond;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::Bond&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::Bond&);

namespace SireCAS
{
class Expression;
}

namespace SireMM
{

/** This class provides a molecule view to a bond */
class SIREMM_EXPORT Bond
    : public SireBase::ConcreteProperty<Bond, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const Bond&);
friend QDataStream& ::operator>>(QDataStream&, Bond&);

public:
    Bond();
    Bond(const SireMol::Atom &atom0, const SireMol::Atom &atom1);
    Bond(const SireMol::MoleculeView &molview,
         const SireMol::AtomID &atom0, const SireMol::AtomID &atom1);
    Bond(const SireMol::MoleculeData &moldata,
         const SireMol::AtomID &atom0, const SireMol::AtomID &atom1);
    Bond(const SireMol::MoleculeData &moldata, const SireMol::BondID &bond);

    Bond(const Bond &other);
    virtual ~Bond();

    static const char* typeName();

    virtual const char* what() const
    {
        return Bond::typeName();
    }

    virtual Bond* clone() const
    {
        return new Bond(*this);
    }

    Bond& operator=(const Bond &bond);

    bool operator==(const Bond &other) const;
    bool operator!=(const Bond &other) const;

    SireMol::MolViewPtr toSelector() const;

    QString toString() const;

    SireMol::Atom atom0() const;
    SireMol::Atom atom1() const;

    SireMol::BondID ID() const;

    bool isEmpty() const;
    bool selectedAll() const;

    SireMol::AtomSelection selection() const;

    bool hasProperty(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key) const;
    bool hasMetadata(const SireBase::PropertyName &key,
                     const SireBase::PropertyName &metakey) const;

    QStringList propertyKeys() const;
    QStringList metadataKeys() const;
    QStringList metadataKeys(const SireBase::PropertyName &key) const;

    SireBase::Properties properties() const;

    SireMol::Evaluator evaluate() const;
    SireMol::Mover<Bond> move() const;

    const Property& property(const SireBase::PropertyName &key) const;
    const Property& property(const SireBase::PropertyName &key,
                             const Property &default_value) const;

    SireUnits::Dimension::Length length() const;
    SireUnits::Dimension::Length length(const SireBase::PropertyMap &map) const;

    SireCAS::Expression potential() const;
    SireCAS::Expression potential(const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

protected:
    bool hasConnectivity() const;
    const SireMol::Connectivity& getConnectivity() const;

    /** The ID of the bond (holding AtomIdx IDs) */
    SireMol::BondID bnd;
};

}

Q_DECLARE_METATYPE( SireMM::Bond )
Q_DECLARE_METATYPE(SireMol::Mover<SireMM::Bond>);

SIRE_EXPOSE_CLASS( SireMM::Bond )
SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::Bond>, SireMol::Mover_Bond_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::Bond>;
}

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
