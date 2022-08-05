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

#ifndef SIREMM_ANGLE_H
#define SIREMM_ANGLE_H

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
class Angle;
}

SIREMM_EXPORT QDataStream& operator<<(QDataStream&, const SireMM::Angle&);
SIREMM_EXPORT QDataStream& operator>>(QDataStream&, SireMM::Angle&);

namespace SireCAS
{
class Expression;
}

namespace SireMM
{

/** This class provides a molecule view to an angle */
class SIREMM_EXPORT Angle
    : public SireBase::ConcreteProperty<Angle, SireMol::MoleculeView>
{

friend QDataStream& ::operator<<(QDataStream&, const Angle&);
friend QDataStream& ::operator>>(QDataStream&, Angle&);

public:
    Angle();
    Angle(const SireMol::Atom &atom0, const SireMol::Atom &atom1,
          const SireMol::Atom &atom2);
    Angle(const SireMol::MoleculeView &molview,
          const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
          const SireMol::AtomID &atom2);
    Angle(const SireMol::MoleculeData &moldata,
          const SireMol::AtomID &atom0, const SireMol::AtomID &atom1,
          const SireMol::AtomID &atom2);
    Angle(const SireMol::MoleculeData &moldata, const SireMol::AngleID &angle);

    Angle(const Angle &other);
    virtual ~Angle();

    static const char* typeName();

    virtual const char* what() const
    {
        return Angle::typeName();
    }

    virtual Angle* clone() const
    {
        return new Angle(*this);
    }

    Angle& operator=(const Angle &other);

    bool operator==(const Angle &other) const;
    bool operator!=(const Angle &other) const;

    SireMol::MolViewPtr toSelector() const;

    QString toString() const;

    SireMol::Atom atom0() const;
    SireMol::Atom atom1() const;
    SireMol::Atom atom2() const;

    SireMol::AngleID ID() const;

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
    SireMol::Mover<Angle> move() const;

    const Property& property(const SireBase::PropertyName &key) const;
    const Property& property(const SireBase::PropertyName &key,
                             const Property &default_value) const;

    SireUnits::Dimension::Angle size() const;
    SireUnits::Dimension::Angle size(const SireBase::PropertyMap &map) const;

    SireCAS::Expression potential() const;
    SireCAS::Expression potential(const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::MolarEnergy energy() const;
    SireUnits::Dimension::MolarEnergy energy(
                            const SireBase::PropertyMap &map) const;

protected:
    bool hasConnectivity() const;
    const SireMol::Connectivity& getConnectivity() const;

    /** The ID of the angle (holding AtomIdx IDs) */
    SireMol::AngleID ang;
};

}

Q_DECLARE_METATYPE( SireMM::Angle )
Q_DECLARE_METATYPE(SireMol::Mover<SireMM::Angle>);

SIRE_EXPOSE_CLASS( SireMM::Angle )
SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMM::Angle>, SireMol::Mover_Angle_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "SireMol/mover.hpp"

namespace SireMol
{
    template class SireMol::Mover<SireMM::Angle>;
}

#endif //SIRE_INSTANTIATE_TEMPLATES

SIRE_END_HEADER

#endif
