/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include "SireMaths/vector.h"

#include "atomproperty.hpp"

#include "atombeads.h"
#include "atomcharges.h"
#include "atomelements.h"
#include "atomenergies.h"
#include "atommasses.h"
#include "atomforces.h"
#include "atomvelocities.h"
#include "atompolarisabilities.h"
#include "atomradii.h"

#include "SireError/errors.h"

using namespace SireMol;
using namespace SireBase;

///////
/////// Implementation of AtomProp
///////

static const RegisterMetaType<AtomProp> r_atomprop(MAGIC_ONLY,
                                                   "SireMol::AtomProp");
                                                   
/** Serialise to a binary datastream */
QDataStream &operator<<(QDataStream &ds, const AtomProp &atomprop)
{
    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, AtomProp &atomprop)
{
    return ds;
}

AtomProp::AtomProp() : MolViewProperty()
{}

AtomProp::AtomProp(const AtomProp &other) : MolViewProperty(other)
{}

AtomProp::~AtomProp()
{}

void AtomProp::throwIncorrectNumberOfAtoms(int nats, int ntotal) const
{
    throw SireError::incompatible_error( QObject::tr(
        "This is incompatible as the number of atoms is %1, but the "
        "number of properties is %2.")
            .arg(ntotal).arg(nats), CODELOC );
}

void AtomProp::throwIncorrectNumberOfSelectedAtoms(int nats, int nselected) const
{
    throw SireError::incompatible_error( QObject::tr(
        "This is incompatible as the number of selected atoms is %1, but the "
        "number of properties is %2.")
            .arg(nselected).arg(nats), CODELOC );
}

////////
//////// Lets explicitly instantiate other AtomProperty types...
////////

namespace SireBase
{
    template class PackedArray2D<SireMol::BeadNum>;
    template class PackedArray2D<SireUnits::Dimension::MolarMass>;
    template class PackedArray2D<SireUnits::Dimension::MolarEnergy>;
    template class PackedArray2D<SireMol::Element>;
    template class PackedArray2D<SireUnits::Dimension::Charge>;
    template class PackedArray2D<SireMol::Velocity3D>;
    template class PackedArray2D<SireMol::Force3D>;
    template class PackedArray2D<SireUnits::Dimension::Volume>;
    template class PackedArray2D<SireUnits::Dimension::Length>;
}

namespace SireMol
{
    template class AtomProperty<QString>;
    template class AtomProperty<qint64>;
    template class AtomProperty<double>;
    template class AtomProperty<QVariant>;
    
    template class AtomProperty<BeadNum>;
    template class AtomProperty<SireUnits::Dimension::MolarMass>;
    template class AtomProperty<SireUnits::Dimension::MolarEnergy>;
    template class AtomProperty<Element>;
    template class AtomProperty<SireUnits::Dimension::Charge>;

    template class AtomProperty<SireMol::Velocity3D>;

    template class AtomProperty<SireMol::Force3D>;
    
    template class AtomProperty<SireUnits::Dimension::Volume>;
    
    template class AtomProperty<SireUnits::Dimension::Length>;
}

static const RegisterMetaType<AtomStringProperty> r_atomstring;
static const RegisterMetaType<AtomIntProperty> r_atomint;
static const RegisterMetaType<AtomFloatProperty> r_atomfloat;
static const RegisterMetaType<AtomVariantProperty> r_atomvariant;
static const RegisterMetaType<BeadNum> r_atombeads;
static const RegisterMetaType<AtomCharges> r_atomcharges;
static const RegisterMetaType<AtomEnergies> r_atomenergies;
static const RegisterMetaType<AtomMasses> r_atommasses;
static const RegisterMetaType<AtomForces> r_atomforces;
static const RegisterMetaType<AtomVelocities> r_atomvelocities;
static const RegisterMetaType<AtomElements> r_atomelements;
static const RegisterMetaType<AtomPolarisabilities> r_atompols;
static const RegisterMetaType<AtomRadii> r_atomradii;
