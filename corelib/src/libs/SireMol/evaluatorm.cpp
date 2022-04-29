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

#include "evaluatorm.h"
#include "atomcoords.h"
#include "atommasses.h"
#include "atomcharges.h"
#include "atomelements.h"
#include "atommatcher.h"
#include "atommatchers.h"
#include "bondid.h"
#include "angleid.h"
#include "dihedralid.h"
#include "connectivity.h"
#include "molecule.h"
#include "mover.hpp"
#include "editor.hpp"
#include "selectormol.h"

#include "SireVol/coordgroup.h"

#include "SireMaths/sphere.h"
#include "SireMaths/axisset.h"
#include "SireMaths/line.h"
#include "SireMaths/triangle.h"
#include "SireMaths/torsion.h"
#include "SireMaths/accumulator.h"

#include "SireBase/errors.h"
#include "SireMol/errors.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

using namespace SireMol;
using namespace SireMaths;
using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<EvaluatorM> r_eval;

SIREMOL_EXPORT QDataStream& operator<<(QDataStream &ds, const EvaluatorM &eval)
{
    writeHeader(ds, r_eval, 1);

    SharedDataStream sds(ds);

    sds << eval.vws << static_cast<const Property&>(eval);

    return ds;
}

SIREMOL_EXPORT QDataStream& operator>>(QDataStream &ds, EvaluatorM &eval)
{
    VersionID v = readHeader(ds, r_eval);

    if (v == 1)
    {
        SharedDataStream sds(ds);
        sds >> eval.vws >> static_cast<Property&>(eval);
    }
    else
        throw version_error(v, "1", r_eval, CODELOC);

    return ds;
}

EvaluatorM::EvaluatorM() : ConcreteProperty<EvaluatorM, Property>()
{}

EvaluatorM::EvaluatorM(const SelectorMol &mols) : ConcreteProperty<EvaluatorM, Property>()
{
    this->vws.reserve(mols.count());

    for (const auto &mol : mols)
    {
        this->vws.append(PartialMolecule(mol));
    }
}

EvaluatorM::EvaluatorM(const EvaluatorM &other)
           : ConcreteProperty<EvaluatorM, Property>(other), vws(other.vws)
{}

EvaluatorM::~EvaluatorM()
{}

const char* EvaluatorM::typeName()
{
    return QMetaType::typeName( qMetaTypeId<EvaluatorM>() );
}

EvaluatorM& EvaluatorM::operator=(const EvaluatorM &other)
{
    if (this != &other)
    {
        vws = other.vws;
        Property::operator=(other);
    }

    return *this;
}

bool EvaluatorM::operator==(const EvaluatorM &other) const
{
    return vws == other.vws;
}

bool EvaluatorM::operator!=(const EvaluatorM &other) const
{
    return not EvaluatorM::operator==(other);
}

int EvaluatorM::nAtoms() const
{
    int n = 0;

    for (const auto &view : this->vws)
    {
        n += view.nAtoms();
    }

    return n;
}

int EvaluatorM::nMolecules() const
{
    return this->vws.count();
}

bool EvaluatorM::isEmpty() const
{
    return this->vws.isEmpty();
}

QString EvaluatorM::toString() const
{
    return QObject::tr("EvaluatorM( nMolecules=%1 nAtoms=%2 )")
            .arg(this->nMolecules())
            .arg(this->nAtoms());
}

SireUnits::Dimension::MolarMass EvaluatorM::mass() const
{
    return this->mass(PropertyMap());
}

SireUnits::Dimension::MolarMass EvaluatorM::mass(const SireBase::PropertyMap &map) const
{
    if (this->isEmpty())
        return SireUnits::Dimension::MolarMass(0);

    auto m = this->vws[0].evaluate().mass(map);

    for (int i=1; i<this->vws.count(); ++i)
    {
        m += this->vws[i].evaluate().mass(map);
    }

    return m;
}

SireUnits::Dimension::Charge EvaluatorM::charge() const
{
    return this->charge(PropertyMap());
}

SireUnits::Dimension::Charge EvaluatorM::charge(const PropertyMap &map) const
{
    if (this->isEmpty())
        return SireUnits::Dimension::Charge(0);

    auto c = this->vws[0].evaluate().charge(map);

    for (int i=1; i<this->vws.count(); ++i)
    {
        c += this->vws[i].evaluate().charge(map);
    }

    return c;
}

Vector EvaluatorM::center() const
{
    return this->center(PropertyMap());
}

Vector EvaluatorM::center(const PropertyMap &map) const
{
    return this->aaBox(map).center();
}

AABox EvaluatorM::aaBox() const
{
    return this->aaBox(PropertyMap());
}

AABox EvaluatorM::aaBox(const PropertyMap &map) const
{
    throw SireError::incomplete_code();
    return AABox();
}

Sphere EvaluatorM::boundingSphere() const
{
    return this->boundingSphere(PropertyMap());
}

Sphere EvaluatorM::boundingSphere(const PropertyMap &map) const
{
    return this->aaBox(map).boundingSphere();
}

Vector EvaluatorM::centroid() const
{
    return this->centroid(PropertyMap());
}

Vector EvaluatorM::centroid(const PropertyMap &map) const
{
    throw SireError::incomplete_code();
    return Vector();
}

Vector EvaluatorM::centerOfGeometry() const
{
    return this->centerOfGeometry(PropertyMap());
}

Vector EvaluatorM::centerOfGeometry(const PropertyMap &map) const
{
    throw SireError::incomplete_code();
    return Vector();
}

Vector EvaluatorM::centerOfMass() const
{
    return this->centerOfMass(PropertyMap());
}

Vector EvaluatorM::centerOfMass(const PropertyMap &map) const
{
    throw SireError::incomplete_code();
    return Vector();
}
